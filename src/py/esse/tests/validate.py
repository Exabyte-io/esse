import os
import unittest

from esse import ESSE
from esse.settings import SCHEMAS_DIR, EXAMPLES_DIR
from esse.utils import parseIncludeReferenceStatements

esse = ESSE()


class TestCase(unittest.TestCase):

    def __init__(self, example_path):
        super(TestCase, self).__init__(methodName="runTest")
        self.example_path = example_path
        self.schema_path = self.example_path.replace(EXAMPLES_DIR, SCHEMAS_DIR)

    def runTest(self):
        schema = parseIncludeReferenceStatements(self.schema_path)
        example = parseIncludeReferenceStatements(self.example_path)
        esse.validate(example, schema)


class TestResult(unittest.TextTestResult):
    def getDescription(self, test):
        return test.example_path.replace("".join((EXAMPLES_DIR, "/")), "")


def validate_examples():
    """
    Validates examples against their schemas.
    """
    suite = unittest.TestSuite()
    for root, dirs, files in os.walk(EXAMPLES_DIR):
        for file_ in files:
            TestCase(example_path=os.path.join(root, file_)).run()
            suite.addTest(TestCase(example_path=os.path.join(root, file_)))
    unittest.TextTestRunner(verbosity=2, resultclass=TestResult).run(suite)


if __name__ == "__main__":
    validate_examples()
