"""
This script builds the schemas and examples modules.
The schemas module contains all the schemas in the SCHEMAS_DIR.
The examples module contains all the examples in the EXAMPLES_DIR.
The script is aware of the include and reference statements and is using the json_include library to resolve them.
The script creates the python modules containing the schemas and examples and avoids
having the dependency on the filesystem during the calls to the ESSE class.
Usage (NOTE: "-e" in pip install is required for the file references to be resolved):
    virtualenv .venv
    source .venv/bin/activate
    pip install -e ."[tests]"
    python build_schemas.py
After that, check for the files `src/py/mat3ra/esse/data/schemas.py` and `src/py/mat3ra/esse/data/examples.py`.
"""
import os
import json
import yaml

from mat3ra.esse.utils import parse_include_reference_statements_by_dir

TOP_DIR = os.path.dirname(__file__)
SCHEMAS = parse_include_reference_statements_by_dir(TOP_DIR)
# Wrap examples into {"data": example, "path": schema_id} format to relate them to schemas.
EXAMPLES = parse_include_reference_statements_by_dir(TOP_DIR, True)

with open(os.path.join(TOP_DIR, "manifest/properties.yaml")) as f:
    PROPERTIES_MANIFEST = yaml.load(f.read(), Loader=yaml.FullLoader)
    RESULTS = [k for k, v in PROPERTIES_MANIFEST.items() if v.get("isResult")]

with open("src/py/mat3ra/esse/data/schemas.py", "w") as f:
    f.write(f"import json; SCHEMAS = json.loads(json.dumps({SCHEMAS}))")

with open("src/py/mat3ra/esse/data/examples.py", "w") as f:
    f.write(f"import json; EXAMPLES = json.loads(json.dumps({EXAMPLES}))")

with open("src/py/mat3ra/esse/data/properties.py", "w") as f:
    content = (f"import json\n" +
               f"PROPERTIES_MANIFEST = json.loads(json.dumps({PROPERTIES_MANIFEST}))\n" +
               f"RESULTS = json.loads(json.dumps({RESULTS}))\n")
    f.write(content)
