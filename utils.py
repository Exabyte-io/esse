import os
import json
import ordereddict

from esse import SCHEMA_DR, EXAMPLE_DR


def dump_json_file(path_, sort_keys=True):
    """
    Dumps a given JSON file. The keys are sorted if `sort_keys` is set, otherwise the order is preserved.

    Args:
        path_ (str): path to JSON file.
        sort_keys (bool): whether to sort keys. Defaults to True.
    """
    with open(path_, "r")as f:
        content = json.loads(f.read(), object_pairs_hook=ordereddict.OrderedDict)
    with open(path_, "w+")as f:
        f.write("".join((json.dumps(content, separators=(',', ': '), indent=4, sort_keys=sort_keys), "\n")))


def dump_examples():
    """
    re-dumps examples:
        - fix space and indentation
        - sort keys
    """
    for root, dirs, files in os.walk(EXAMPLE_DR):
        for file_ in files:
            dump_json_file(os.path.join(root, file_))


def dump_schemas():
    """
    re-dumps schemas:
        - fix space and indentation
    """
    for root, dirs, files in os.walk(SCHEMA_DR):
        for file_ in files:
            dump_json_file(os.path.join(root, file_), False)


if __name__ == "__main__":
    dump_schemas()
    dump_examples()
