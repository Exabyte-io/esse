/* eslint-disable class-methods-use-this */
import Ajv from "ajv";
import fs from "fs";
import yaml from "js-yaml";
import mergeAllOf from "json-schema-merge-allof";
import path from "path";

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
        const ajv = new Ajv({ allErrors: true, allowUnionTypes: true, discriminator: true });
        const isValid = ajv.validate(schema, example);

        if (!isValid) {
            console.error(ajv.errors);
        }

        return isValid;
    };

    writeResolvedSchemas(subfolder, skipMergeAllOff = false) {
        this.schemas.forEach((s) => {
            let mergedSchema = s;
            if (!skipMergeAllOff) {
                mergedSchema = mergeAllOf(s, {
                    resolvers: { defaultResolver: mergeAllOf.options.resolvers.title },
                });
            }
            const id_as_path = mergedSchema.$id.replace("-", "_");
            const full_path = `${subfolder}/schema/${id_as_path}.json`;
            fs.mkdirSync(path.dirname(full_path), { recursive: true });
            fs.writeFileSync(full_path, JSON.stringify(mergedSchema, null, 4), "utf8");
        });
    }

    writeResolvedExamples(subfolder) {
        this.wrappedExamples.forEach((e) => {
            const id_as_path = e.path.replace("-", "_");
            const full_path = `${subfolder}/example/${id_as_path}.json`;
            fs.mkdirSync(path.dirname(full_path), { recursive: true });
            fs.writeFileSync(full_path, JSON.stringify(e.data, null, 4), "utf8");
        });
    }
}
