import os
import json
import json_include


def parse_include_reference_statements(file_path):
    """
    Resolves `include` and `$ref` statements.

    Args:
        file_path (str): file to parse.

    Returns:
         dict|list
    """
    dir_name = os.path.dirname(file_path)
    base_name = os.path.basename(file_path)
    return json.loads(json_include.build_json(dir_name, base_name))


def parse_include_reference_statements_by_dir(top_dir, is_examples=False):
    """
    Resolves `include` and `$ref` statements for all the JSON files inside a given directory.

    Args:
        dir_path (str): directory to parse.

    Returns:
         dict|list
         :param top_dir: Top-level directory containing "schema" and "example" directories.
         :param is_examples: Whether to wrap and add the path to the returned data.
    """
    data = []
    dir_path = os.path.join(top_dir, "schema") if not is_examples else os.path.join(top_dir, "example")
    for root, dirs, files in os.walk(dir_path):
        for file_ in files:
            if os.path.splitext(file_)[1] == ".json":
                file_path = os.path.join(root, file_)
                if is_examples:
                    schema_id = (
                        file_path.replace("example", "schema")
                        .replace(top_dir, "")
                        .replace("schema", "")
                        .replace(".json", "")
                        .strip("/")
                    )
                    config = {"data": parse_include_reference_statements(file_path), "path": schema_id}
                else:
                    config = parse_include_reference_statements(file_path)
                data.append(config)
    return data
