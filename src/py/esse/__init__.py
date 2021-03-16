import jsonschema

from esse.utils import parseIncludeReferenceStatementsByDir
from esse.settings import SCHEMAS_DIR, EXAMPLES_DIR, PROPERTIES_MANIFEST
# from esse.functionals import generate_dft_unit_functionals

## Note: We're disabling this for now to avoid a bug discovered in deployment - files are not necessarily writeable
## in the installed context, if the calling code is executed with a different user than the one that owns the install
## directory.  The function can be run manually with the `generate_dft_unit_functionals` command line script generated
## in setup.py when esse is installed.
# generate_dft_unit_functionals()
SCHEMAS = parseIncludeReferenceStatementsByDir(SCHEMAS_DIR)
EXAMPLES = parseIncludeReferenceStatementsByDir(EXAMPLES_DIR)


class ESSE(object):
    """
    Exabyte Source of Schemas and Examples class.
    """

    def __init__(self):
        self.schemas = SCHEMAS
        self.examples = EXAMPLES

    def get_schema_by_id(self, schemaId):
        return next((s for s in SCHEMAS if s.get("schemaId") == schemaId), None)

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

    def get_property_manifest(self, property_):
        """
        Returns the manifest for a given property.

        Args:
            property_ (str): property name.

        Returns:
             dict
        """
        return PROPERTIES_MANIFEST.get(property_, {})
