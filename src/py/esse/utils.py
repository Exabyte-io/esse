import os
import json
import json_include

from slugify import slugify
from collections import OrderedDict

from esse.settings import SCHEMAS_DIR, EXAMPLES_DIR


def parseIncludeReferenceStatements(file_path):
    """
    Resolves `include` and `$ref` statements.

    Args:
        file_path (str): file to parse.

    Returns:
         dict|list
    """
    dirName = os.path.dirname(file_path)
    baseName = os.path.basename(file_path)
    return json.loads(json_include.build_json(dirName, baseName))


def parseIncludeReferenceStatementsByDir(dir_path):
    """
    Resolves `include` and `$ref` statements for all the JSON files inside a given directory.

    Args:
        dir_path (str): directory to parse.

    Returns:
         dict|list
    """
    data = []
    for root, dirs, files in os.walk(dir_path):
        for file_ in files:
            file_path = os.path.join(root, file_)
            data.append(parseIncludeReferenceStatements(file_path))
    return data


def read_json_file(path_):
    """
    Reads the given JSON file into an ordered dictionary.

    Args:
        path_ (str): path to JSON file.

    Returns:
        OrderedDict
    """
    with open(path_, "r") as f:
        return json.loads(f.read(), object_pairs_hook=OrderedDict)


def dump_json_file(path_, content, sort_keys=True):
    """
    Dumps a given JSON content. The keys are sorted if `sort_keys` is set.

    Args:
        path_ (str): path to JSON file.
        content (dict): JSON file content.
        sort_keys (bool): whether to sort keys. Defaults to True.
    """
    with open(path_, "w+") as f:
        f.write("".join((json.dumps(content, separators=(',', ': '), indent=4, sort_keys=sort_keys), "\n")))


def refactor_examples():
    """
    Refactors examples:
        - fix space and indentation
        - sort keys
    """
    for root, dirs, files in os.walk(EXAMPLES_DIR):
        for file_ in files:
            path_ = os.path.join(root, file_)
            content = read_json_file(path_)
            dump_json_file(path_, content, True)


def refactor_schemas():
    """
    Refactors schemas:
        - fix space and indentation
        - preserve keys order
    """
    for root, dirs, files in os.walk(SCHEMAS_DIR):
        for file_ in files:
            set_schema_id(os.path.join(root, file_))


def set_schema_id(path_):
    """
    Slugifies the schema path and uses it as schema ID.

    Args:
        path_ (str): path to JSON file.
    """
    content = read_json_file(path_)
    if not content.get("$schema"): return  # do not add ID to non-schema files
    if content.get("schemaId"): del content["schemaId"]
    schema_id = slugify(path_.replace("{}/".format(SCHEMAS_DIR), '').replace(".json", ""))
    content = OrderedDict(list(OrderedDict({"schemaId": schema_id}).items()) + list(content.items()))
    dump_json_file(path_, content, False)
