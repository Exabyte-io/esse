import os
import json

from slugify import slugify
from collections import OrderedDict

from esse import SCHEMA_DR, EXAMPLE_DR


def read_json_file(path_):
    """
    Reads the given JSON file into an ordered dictionary.

    Args:
        path_ (str): path to JSON file.

    Returns:
        OrderedDict
    """
    with open(path_, "r")as f:
        return json.loads(f.read(), object_pairs_hook=OrderedDict)


def dump_json_file(path_, content, sort_keys=True):
    """
    Dumps a given JSON content. The keys are sorted if `sort_keys` is set.

    Args:
        path_ (str): path to JSON file.
        content (dict): JSON file content.
        sort_keys (bool): whether to sort keys. Defaults to True.
    """
    with open(path_, "w+")as f:
        f.write("".join((json.dumps(content, separators=(',', ': '), indent=4, sort_keys=sort_keys), "\n")))


def dump_examples():
    """
    Dumps examples:
        - fix space and indentation
        - sort keys
    """
    for root, dirs, files in os.walk(EXAMPLE_DR):
        for file_ in files:
            path_ = os.path.join(root, file_)
            content = read_json_file(path_)
            dump_json_file(path_, content, True)


def dump_schemas():
    """
    Dumps schemas:
        - fix space and indentation
        - preserve keys order
    """
    for root, dirs, files in os.walk(SCHEMA_DR):
        for file_ in files:
            set_schema_id(os.path.join(root, file_))


def set_schema_id(path_):
    """
    Slugifies the schema path and uses it as schema ID.

    Args:
        path_ (str): path to JSON file.
    """
    content = read_json_file(path_)
    if content.get("$id"): del content["$id"]
    schema_id = slugify(path_.replace("{}/".format(SCHEMA_DR), '').replace(".json", ""))
    content = OrderedDict(list(OrderedDict({"$id": schema_id}).items()) + list(content.items()))
    dump_json_file(path_, content, False)


if __name__ == "__main__":
    dump_schemas()
    dump_examples()
