import os
import yaml

DIR = os.path.dirname(__file__)
SCHEMAS_DIR = os.path.join(DIR, "data/schema")
EXAMPLES_DIR = os.path.join(DIR, "data/example")
MANIFEST_DR = os.path.join(DIR, "data/manifest")

with open(os.path.join(MANIFEST_DR, "properties.yaml")) as f:
    PROPERTIES_MANIFEST = yaml.load(f.read(), Loader=yaml.FullLoader)
    RESULTS = [k for k, v in PROPERTIES_MANIFEST.items() if v.get("isResult")]
