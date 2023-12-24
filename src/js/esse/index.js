import Ajv from "ajv";

import { buildSchemaDefinitions } from "./schemaUtils";
import { EXAMPLES_DIR, SCHEMAS_DIR } from "./settings";
import { parseIncludeReferenceStatementsByDir } from "./utils";

const SCHEMAS = parseIncludeReferenceStatementsByDir(SCHEMAS_DIR);
const EXAMPLES = parseIncludeReferenceStatementsByDir(EXAMPLES_DIR, true);

export class ESSE {
    constructor() {
        this.schemas = SCHEMAS;
        this.wrappedExamples = EXAMPLES;
        this.examples = EXAMPLES.map((example) => example.data);
    }

    getSchemaById(schemaId) {
        return this.schemas.find((schema) => schema.$id === schemaId);
    }

    /**
     * Validates a given example against the schema.
     * @param example {Object|Array} example to validate.
     * @param schema {Object} schema to validate the example with.
     * @returns {boolean} whether example is valid.
     */
    validate = (example, schema) => {
        const ajv = new Ajv({ allErrors: true });
        const isValid = ajv.validate(schema, example);

        if (!isValid) {
            console.error(ajv.errors);
        }

        return isValid;
    };

    buildGlobalSchema() {
        return {
            $id: "esse-global-schema",
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "Global schema",
            type: "object",
            definitions: buildSchemaDefinitions(this.schemas),
        };
    }
}
