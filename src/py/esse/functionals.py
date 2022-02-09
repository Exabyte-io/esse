import os
import sys
import json
import errno
from string import Template

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
REL_DIR = "schema/models_directory/pb/qm/dft"
UNIT_FUNCTIONALS_PATH = os.path.join(BASE_DIR, REL_DIR)
PROTOTYPE_FILENAME = os.path.join(UNIT_FUNCTIONALS_PATH, "dft_unit_functionals_proto.json")
UNIT_FILENAME = os.path.join(UNIT_FUNCTIONALS_PATH, "dft_unit_functionals.json")
DATA_FILENAME = os.path.join(BASE_DIR, "src/py/esse/data", REL_DIR, "dft_unit_functionals.json")

SCHEMA_WITH_PROPERTIES_TEMPLATE = Template("""{
"properties": {
    "name": {
        "enum": [
            "$name"
        ]
    },
    "slug": {
        "enum": [
            "$slug"
        ]
    },
    "type": {
        "enum": [
            "$type"
        ]
    }
}
}""")


def remove(fpath):
    try:
        os.unlink(fpath)
    except FileNotFoundError:
        pass
    try:
        os.remove(fpath)
    except FileNotFoundError:
        pass


def generate_dft_unit_functionals():
    """
    Generate list of functionals suitable for validation by 'oneOf'.

    The targeted JSON file has a convoluted form which makes it inconvenient
    to read and prone to errors. Thus, it is generated automatically from
    a prototype file that is easier to maintain.
    """
    with open(PROTOTYPE_FILENAME, "r") as f:
        proto = json.loads(f.read())

    del proto["description"]
    output = {}
    for rung, configs in proto.items():
        o = {"oneOf": []}
        for config in configs:
            o["oneOf"].append(json.loads(SCHEMA_WITH_PROPERTIES_TEMPLATE.substitute(config)))
        output[rung] = o

    remove(UNIT_FILENAME)
    remove(DATA_FILENAME)

    with open(UNIT_FILENAME, "w") as f:
        f.write("".join((json.dumps(output, separators=(',', ': '), indent=4, sort_keys=True), "\n")))


if __name__ == "__main__":
    # add comment
    generate_dft_unit_functionals()
