import Ajv from "ajv";
import keyBy from "lodash/keyBy";

import { EXAMPLES_DIR, SCHEMAS_DIR } from "./settings";
import {
    makeSchemaId,
    makeSchemaRef,
    mapObjectDeep,
    parseIncludeReferenceStatementsByDir,
} from "./utils";

const SCHEMAS = parseIncludeReferenceStatementsByDir(SCHEMAS_DIR);
const EXAMPLES = parseIncludeReferenceStatementsByDir(EXAMPLES_DIR);

export class ESSE {
    constructor() {
        this.schemas = SCHEMAS;
        this.examples = EXAMPLES;
    }

    getSchemaById(schemaId) {
        return this.schemas.find((schema) => schema.schemaId === schemaId);
    }

    /**
     * Validates a given example against the schema.
     * @param example {Object|Array} example to validate.
     * @param schema {Object} schema to validate the example with.
     * @returns {boolean} whether example is valid.
     */
    validate = (example, schema) => {
        const ajv = new Ajv({ allErrors: true });
        return ajv.validate(schema, example);
    };

    buildSchemaDefinitions() {
        const schemas = this.schemas.map((schema) => {
            return mapObjectDeep(schema, (value) => {
                if (typeof value === "object" && value.schemaId) {
                    return makeSchemaRef(value.schemaId);
                }
            });
        });

        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "Global schema",
            // allOf: this.schemas.map(({ schemaId }) => makeRef(schemaId)),
            type: "object",
            definitions: keyBy(schemas, ({ schemaId }) => makeSchemaId(schemaId)),
        };
    }
}
