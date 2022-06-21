#!/usr/bin/env python

import json
from pathlib import Path


ESSE_SCHEMA_DIR = "schema"


def has_schema_id(filepath: Path, identifier="schemaId") -> bool:
    with filepath.open() as f:
        schema = json.load(f)
    return identifier in schema


def get_schema_dict(path: Path) -> dict:
    with path.open(mode="r") as f:
        schema = json.load(f)
    return schema


def write_schema(path: Path, schema: dict) -> None:
    with path.open(mode="w") as out:
        json.dump(schema, out, indent=4)
        out.write("\n")


def schema_id_from_path(path: Path) -> str:
    rel_path = path.relative_to(ESSE_SCHEMA_DIR)
    return (rel_path.parent / path.stem).as_posix()


def main() -> None:
    base_path = Path(ESSE_SCHEMA_DIR)
    schema_paths = [s for s in base_path.rglob("*.json") if has_schema_id(s)]
    print(f"Found {len(schema_paths)} schema IDs.")

    for i, path in enumerate(schema_paths):
        schema = get_schema_dict(path)
        schema["schemaId"] = schema_id_from_path(path)
        write_schema(path, schema)
    print("Schema ID replacement complete.")


if __name__ == "__main__":
    main()
