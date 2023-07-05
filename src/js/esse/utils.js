import file from "file";
import deref from "json-schema-deref-sync";
import path from "path";

import { JSONInclude } from "../json_include";

/**
 * Resolves `include` and `$ref` statements.
 * @param filePath {String} file to parse.
 */
export function parseIncludeReferenceStatements(filePath) {
    const jsonResolver = new JSONInclude();
    const parsed = jsonResolver.parseIncludeStatements(filePath);
    const dirPath = path.dirname(filePath);
    let dereferenced = deref(parsed, { baseFolder: dirPath });
    // handle circular references and use non-dereferenced source
    if (dereferenced instanceof Error && dereferenced.message === "Circular self reference") {
        dereferenced = parsed;
    }
    return dereferenced;
}

/**
 * Resolves `include` and `$ref` statements for all the JSON files inside a given directory.
 * @param dirPath {String} directory to parse.
 */
export function parseIncludeReferenceStatementsByDir(dirPath) {
    const data = [];
    file.walkSync(dirPath, (dirPath_, dirs_, files_) => {
        files_.forEach((file_) => {
            const filePath = path.join(dirPath_, file_);
            data.push(parseIncludeReferenceStatements(filePath));
        });
    });
    return data;
}

export function mapObjectDeep(object, mapValue) {
    if (typeof object !== "object") {
        return object;
    }

    if (Array.isArray(object)) {
        return object.map(
            (innerValue) => mapValue(innerValue) || mapObjectDeep(innerValue, mapValue),
        );
    }

    const entries = Object.entries(object).map(([key, value]) => {
        return [key, mapValue(value) || mapObjectDeep(value, mapValue)];
    });

    return Object.fromEntries(entries);
}

export function makeSchemaId(schemaId) {
    return schemaId.replace(/\//g, "-");
}

export function makeSchemaRef(schemaId) {
    return { $ref: `#/definitions/${makeSchemaId(schemaId)}` };
}
