import json
from string import Template
import os
from esse.utils import read_json_file, dump_json_file

DIR = os.path.dirname(__file__)
UNIT_FUNCTIONALS_PATH = os.path.join(DIR, "data/schema/models_directory/pb/qm/dft")
PROTOTYPE_FILENAME = os.path.join(UNIT_FUNCTIONALS_PATH, "dft_unit_functionals_proto.json")
UNIT_FILENAME = os.path.join(UNIT_FUNCTIONALS_PATH, "dft_unit_functionals.json")

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

def generate_dft_unit_functionals():
    """
    Generate list of functionals suitable for validation by 'oneOf'.

    The targeted JSON file has a convoluted form which makes it inconvenient
    to read and prone to errors. Thus, it is generated automatically from
    a prototype file that is easier to maintain.
    """
    proto = read_json_file(PROTOTYPE_FILENAME)
    del proto["description"]
    output = {}
    for rung, configs in proto.items():
        o = {"oneOf": []}
        for config in configs:
            o["oneOf"].append(json.loads(SCHEMA_WITH_PROPERTIES_TEMPLATE.substitute(config)))
        output[rung] = o

    dump_json_file(UNIT_FILENAME, output)
