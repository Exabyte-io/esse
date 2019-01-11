import os
import json
import jsonschema
import json_include

from settings import SCHEMAS_DR, EXAMPLES_DIR, PROPERTIES_MANIFEST


class ESSE(object):
    """
    Exabyte Source of Schemas and Examples class.
    """

    def __init__(self):
        self.schemas = self.parseIncludeReferenceStatementsByDir(SCHEMAS_DR)
        self.examples = self.parseIncludeReferenceStatementsByDir(EXAMPLES_DIR)

    def get_schema_by_id(self, schemaId):
        return next((s for s in self.schemas if s["schemaId"] == schemaId), None)

    def validate(self, example, schema):
        """
        Validates a given example against the schema.

        Args:
            example (dict|list): example to validate.
            schema (dict): schema to validate the example with.

        Raises:
            jsonschema.exceptions.ValidationError
        """
        jsonschema.validate(example, schema)

    def parseIncludeReferenceStatements(self, file_path):
        """
        Resolves `include` and `$ref` statements.

        Args:
            file_path (str): file to parse.

        Returns:
             dict|list
        """
        dirName = os.path.dirname(file_path)
        baseName = os.path.basename(file_path)
        return json.loads(json_include.build_json(dirName, baseName))

    def parseIncludeReferenceStatementsByDir(self, dir_path):
        """
        Resolves `include` and `$ref` statements for all the JSON files inside a given directory.

        Args:
            dir_path (str): directory to parse.

        Returns:
             dict|list
        """
        data = []
        for root, dirs, files in os.walk(dir_path, followlinks=True):
            for file_ in files:
                file_path = os.path.join(root, file_)
                data.append(self.parseIncludeReferenceStatements(file_path))
        return data

    def get_property_manifest(self, property_):
        """
        Returns the manifest for a given property.

        Args:
            property_ (str): property name.

        Returns:
             dict
        """
        return PROPERTIES_MANIFEST.get(property_, {}).get("defaults", {})
