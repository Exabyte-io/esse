import Ajv from "ajv";
import file from "file";
import path from "path";
import deref from "json-schema-deref-sync";

import {JSONInclude} from "../json_include/index";
import {EXAMPLES_DIR, SCHEMAS_DIR} from "./settings";


export class ESSE {

    constructor() {
        this.jsonResolver = new JSONInclude();
        this.schemas = this.parseIncludeReferenceStatementsByDir(SCHEMAS_DIR);
        this.examples = this.parseIncludeReferenceStatementsByDir(EXAMPLES_DIR);
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
        return ajv.validate(schema, example);
    }

    /**
     * Resolves `include` and `$ref` statements.
     * @param filePath {String} file to parse.
     */
    parseIncludeReferenceStatements(filePath) {
        const parsed = this.jsonResolver.parseIncludeStatements(filePath);
        const dirPath = path.dirname(filePath);
        let dereferenced = deref(parsed, {baseFolder: dirPath});
        // handle circular references and use non-dereferenced source
        if ((dereferenced instanceof Error) && (dereferenced.message === "Circular self reference")) {
            dereferenced = parsed;
        }
        return dereferenced;
    }

    /**
     * Resolves `include` and `$ref` statements for all the JSON files inside a given directory.
     * @param dirPath {String} directory to parse.
     */
    parseIncludeReferenceStatementsByDir(dirPath) {
        const data = [];
        file.walkSync(dirPath, (dirPath_, dirs_, files_) => {
            files_.forEach(file_ => {
                const filePath = path.join(dirPath_, file_);
                data.push(this.parseIncludeReferenceStatements(filePath));
            });
        });
        return data;
    }
}
