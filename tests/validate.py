import os
import json

from esse import ESSE, SCHEMA_DR, EXAMPLE_DR


def sort_json_fields(path_):
    """
    Sorts JSON fields by key.

    Args:
        path_ (str): path to JSON file.
    """
    with open(path_, "r")as f:
        content = json.loads(f.read())
    with open(path_, "w+")as f:
        f.write("".join((json.dumps(content, indent=4, sort_keys=True), "\n")))


def sort_examples():
    """
    Sorts examples.
    """
    for root, dirs, files in os.walk(EXAMPLE_DR):
        for file_ in files:
            sort_json_fields(os.path.join(root, file_))


def sort_schemas():
    """
    Sorts schemas.
    """
    for root, dirs, files in os.walk(SCHEMA_DR):
        for file_ in files:
            sort_json_fields(os.path.join(root, file_))


def validate_examples():
    """
    Validates examples agains their schemas.
    """
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
    sort_examples()
    validate_examples()
