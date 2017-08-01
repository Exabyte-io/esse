import os
import yaml
import jsonschema
import jsonInclude

ESSE_ROOT_DIR = os.path.dirname(__file__)
SCHEMA_DR = os.path.join(ESSE_ROOT_DIR, "schema")
EXAMPLE_DR = os.path.join(ESSE_ROOT_DIR, "example")

with open(os.path.join(ESSE_ROOT_DIR, "manifest", "schemas.yaml")) as f:
    SCHEMAS_MANIFEST = yaml.load(f.read())


class ESSE(object):
    """
    Example and Schema Sources for Exabyte class.
    """

    def get_example(self, schemaId):
        """
        Returns an example for a given schema ID.

        Returns:
             dict
        """
        return self._get_json(os.path.join(ESSE_ROOT_DIR, SCHEMAS_MANIFEST[schemaId]["examplePath"]))

    def get_schema(self, schemaId):
        """
        Returns a schema for a given schema ID.

        Returns:
             dict
        """
        return self._get_json(os.path.join(ESSE_ROOT_DIR, SCHEMAS_MANIFEST[schemaId]["schemaPath"]))

    def get_schema_default_values(self, schemaId):
        """
        Returns default values for a given schema ID.

        Returns:
             dict
        """
        return SCHEMAS_MANIFEST[schemaId]["defaults"]

    def get_schema_ids(self):
        """
        Returns a list of schema IDs.

        Returns:
             list
        """
        return SCHEMAS_MANIFEST.keys()

    def validate(self, schema, instance):
        """
        Validates a given data against the schema.

        Args:
            schema (dict): schema to validate istance against.
            instance (dict): instance to validate.

        Raises:
            jsonschema.exceptions.ValidationError
        """
        jsonschema.Draft4Validator(schema).validate(instance)

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
        return jsonInclude.build_json(dirName, baseName)
