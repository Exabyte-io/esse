import os
import yaml

DIR = os.path.dirname(__file__)
SCHEMA_DR = os.path.join(DIR, "data/schema")
EXAMPLE_DR = os.path.join(DIR, "data/example")
MANIFEST_DR = os.path.join(DIR, "data/manifest")

with open(os.path.join(MANIFEST_DR, "properties.yaml")) as f:
    PROPERTIES_MANIFEST = yaml.load(f.read())
