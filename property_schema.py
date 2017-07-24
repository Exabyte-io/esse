import os
import json
import yaml
import jsonschema
import json_include

THIS_FILE_DIR = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(THIS_FILE_DIR, "manifest", "properties.yaml")) as f:
    PROPERTIES_MANIFEST = yaml.load(f.read())


class PropertySchema(object):
    """
    Property schema class.

    Args:
        name (str): property name.
    """

    def __init__(self, name):
        self._name = name
        json_path = os.path.join(THIS_FILE_DIR, PROPERTIES_MANIFEST[self.name]["schema_path"])
        dirpath = os.path.dirname(json_path)
        filename = os.path.basename(json_path)
        self._schema = json.loads(json_include.build_json_include(dirpath, filename))
        self._defaults = PROPERTIES_MANIFEST[self.name].get('defaults')

    @property
    def schema(self):
        return self._schema

    @property
    def name(self):
        return self._name

    @property
    def defaults(self):
        return self._defaults

    def validate(self, data):
        """
        Validates a given data against the schema.

        Args:
            data (dict): data to validate.

        Raises:
            jsonschema.exceptions.ValidationError
        """
        jsonschema.Draft4Validator(self.schema).validate(data)
