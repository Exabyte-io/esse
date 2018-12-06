import os
import json

from esse import ESSE, SCHEMA_DR, EXAMPLE_DR


def sort_example_fields():
    """
    Sort example fields by key.
    """
    for root, dirs, files in os.walk(EXAMPLE_DR):
        for file_ in files:
            example_path = os.path.join(root, file_)
            with open(example_path, "r")as f:
                content = json.loads(f.read())
            with open(example_path, "w+")as f:
                f.write("".join((json.dumps(content, indent=4, sort_keys=True), "\n")))


def validate_examples():
    es = ESSE()
    for root, dirs, files in os.walk(EXAMPLE_DR):
        for file_ in files:
            example_path = os.path.join(root, file_)
            schema_path = example_path.replace(EXAMPLE_DR, SCHEMA_DR)
            try:
                schema = es._get_json(schema_path)
                example = es._get_json(example_path)
                es.validate(example, schema)
            except:
                print("{} is invalid!".format(example_path))
                raise


if __name__ == "__main__":
    sort_example_fields()
    validate_examples()
