import os
import yaml
import json
import jsonschema
import json_include

ESSE_ROOT_DIR = os.path.dirname(__file__)
SCHEMA_DR = os.path.join(ESSE_ROOT_DIR, "schema")
EXAMPLE_DR = os.path.join(ESSE_ROOT_DIR, "example")

with open(os.path.join(ESSE_ROOT_DIR, "manifest", "properties.yaml")) as f:
    PROPERTIES_MANIFEST = yaml.load(f.read())


class ESSE(object):
    """
    Exabyte Source of Schemas and Examples class.
    """

    def get_schema(self, schemaId):
        return self._get_json(self._get_schema_path(schemaId))

    def validate(self, instance, schema):
        """
        Validates a given data against the schema.

        Args:
            instance (dict): instance to validate.
            schema (dict): schema to validate istance against.

        Raises:
            jsonschema.exceptions.ValidationError
        """
        jsonschema.validate(instance, schema)

    def _get_json(self, path):
        """
        Returns a json with inclusion references resolved.

        Args:
            path (str): path to the json file.

        Returns:
             dict
        """
        dirName = os.path.dirname(path)
        baseName = os.path.basename(path)
        return json.loads(json_include.build_json(dirName, baseName))

    def _get_schema_path(self, schemaId):
        path = PROPERTIES_MANIFEST.get(schemaId, {}).get("path")
        return os.path.join(SCHEMA_DR, path) if path else self._find_file(schemaId, SCHEMA_DR)

    def _find_file(self, name, path):
        for root, dirs, files in os.walk(path, followlinks=True):
            for file_ in files:
                if name in file_:
                    return os.path.join(root, file_)

    def get_property_default_values(self, property_):
        """
        Returns default values for a given property.

        Args:
            property_ (str): property name.

        Returns:
             dict
        """
        return PROPERTIES_MANIFEST.get(property_, {}).get("defaults", {})
