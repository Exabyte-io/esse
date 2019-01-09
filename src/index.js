import Ajv from "ajv";
import prettyjson from "prettyjson";

import {EXAMPLES_DIR, SCHEMAS_DIR} from "./settings";
import {JSONSchemaResolver} from "./resolver/resolver";


export class ESSE {

    constructor() {
        const jsonResolver = new JSONSchemaResolver();
        this.schemas = jsonResolver.resolveDir(SCHEMAS_DIR);
        this.examples = jsonResolver.resolveDir(EXAMPLES_DIR);
    }

    getSchemaById(schemaId) {
        return this.schemas.find(schema => schema.schemaId === schemaId)
    }

    /**
     * Validates a given example against the schema.
     * @param example {Object|Array} example to validate.
     * @param schema {Object} schema to validate the example with.
     * @param printErrors {boolean} whether to print errors.
     * @returns {boolean} whether example is valid.
     */
    validate(example, schema, printErrors = false) {
        const ajv = new Ajv({allErrors: true});
        const valid = ajv.validate(schema, example);
        if (!valid && printErrors) console.log(prettyjson.render(ajv.errors));
        return valid;
    }
}
