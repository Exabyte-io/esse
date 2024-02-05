import Ajv from "ajv";
import fs from "fs";
import yaml from "js-yaml";

import { EXAMPLES_DIR, PROPERTIES_MANIFEST_PATH, SCHEMAS_DIR } from "./settings";
import { parseIncludeReferenceStatementsByDir } from "./utils";

const SCHEMAS = parseIncludeReferenceStatementsByDir(SCHEMAS_DIR);
const EXAMPLES = parseIncludeReferenceStatementsByDir(EXAMPLES_DIR, true);
const PROPERTIES_MANIFEST = yaml.load(
    fs.readFileSync(PROPERTIES_MANIFEST_PATH, { encoding: "utf-8" }),
);
const RESULTS = Object.entries(PROPERTIES_MANIFEST)
    .map((k) => (k[1].isResult ? k[0] : null))
    .filter((x) => x);

export class ESSE {
    constructor(config = {}) {
        this.schemas = config.schemas || SCHEMAS;
        this.wrappedExamples = config.wrappedExamples || EXAMPLES;
        this.examples = this.wrappedExamples.map((example) => example.data);
        this.propertiesManifest = config.propertiesManifest || PROPERTIES_MANIFEST;
        this.results = config.results || RESULTS;
    }

    /**
     * Validates a given example against the schema.
     * @param example {Object|Array} example to validate.
     * @param schema {Object} schema to validate the example with.
     * @returns {boolean} whether example is valid.
     */
    validate = (example, schema) => {
        const ajv = new Ajv({ allErrors: true, allowUnionTypes: true });
        const isValid = ajv.validate(schema, example);

        if (!isValid) {
            console.error(ajv.errors);
        }

        return isValid;
    };
}
