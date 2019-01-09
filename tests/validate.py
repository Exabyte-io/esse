import os
import unittest

from esse import ESSE, SCHEMA_DR, EXAMPLE_DR


class TestCase(unittest.TestCase):

    def __init__(self, example_path):
        super(TestCase, self).__init__(methodName="runTest")
        self.example_path = example_path
        self.schema_path = self.example_path.replace(EXAMPLE_DR, SCHEMA_DR)
        self.esse = ESSE()

    def runTest(self):
        schema = self.esse._get_json(self.schema_path)
        example = self.esse._get_json(self.example_path)
        self.esse.validate(example, schema)


class TestResult(unittest.TextTestResult):
    def getDescription(self, test):
        return test.example_path.replace("".join((EXAMPLE_DR, "/")), "")


def validate_examples():
    """
    Validates examples against their schemas.
    """
    suite = unittest.TestSuite()
    for root, dirs, files in os.walk(EXAMPLE_DR):
        for file_ in files:
            TestCase(example_path=os.path.join(root, file_)).run()
            suite.addTest(TestCase(example_path=os.path.join(root, file_)))
    unittest.TextTestRunner(verbosity=2, resultclass=TestResult).run(suite)


if __name__ == "__main__":
    validate_examples()
