import json
from string import Template
import os
from esse.utils import read_json_file, dump_json_file

DIR = os.path.dirname(__file__)
UNIT_PATH = os.path.join(DIR, "data/schema/models_directory/pb/qm/dft")
PROTO_FILE = os.path.join(UNIT_PATH, "dft_unit_functionals_proto.json")
UNIT_FILE = os.path.join(UNIT_PATH, "dft_unit_functionals.json")

PROP_TEMPLATE = Template("""{
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
    proto = read_json_file(PROTO_FILE)
    del proto["description"]
    output = {}
    for rung, configs in proto.iteritems():
        o = {"oneOf": []}
        for config in configs:
            tmp_string = PROP_TEMPLATE.substitute(config)
            o["oneOf"].append(json.loads(tmp_string))
        output[rung] = o

    dump_json_file(UNIT_FILE, output)
