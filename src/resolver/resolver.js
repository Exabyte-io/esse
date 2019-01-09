import file from "file";
import path from "path";
import fs from "fs-extra";
import deref from "json-schema-deref-sync";

import {safeParseJSON, isInstanceOf} from "./utils";
import {INCLUDE_KEY, INCLUDE_VALUE_PATTERN, OBJECT_ONLY} from "./settings";


export class JSONSchemaResolver {

    constructor() {
        this.JSON_INCLUDE_CACHE = {};
    }

    _getIncludeFileName(value) {
        if ((isInstanceOf(value, "String")) && (value.search(INCLUDE_VALUE_PATTERN) !== -1)) {
            return value.match(INCLUDE_VALUE_PATTERN)[1]
        }
    }

    /**
     * Walks a nested object to resolve include statements.
     * @param obj {Object} Object to traverse.
     * @param dirpath {String} directory from which `obj` is obtained. Include statements are relative to `obj` path.
     */
    _walkObjectToInclude(obj, dirpath) {
        if (isInstanceOf(obj, "Object")) {
            let isIncludeExp = false;
            if (INCLUDE_KEY in obj) {
                const includeName = this._getIncludeFileName(obj[INCLUDE_KEY]);
                if (includeName) {
                    isIncludeExp = true;
                    delete obj[INCLUDE_KEY];
                    if (!(includeName in this.JSON_INCLUDE_CACHE)) {
                        const filePath = path.join(dirpath, includeName);
                        this.JSON_INCLUDE_CACHE[includeName] = this._parseIncludeStatements(filePath);
                    }
                    for (let attr in this.JSON_INCLUDE_CACHE[includeName]) {
                        obj[attr] = this.JSON_INCLUDE_CACHE[includeName][attr];
                    }
                }
            }
            if (isIncludeExp) return;
            for (let key in obj) {
                if (isInstanceOf(obj[key], "Object") || isInstanceOf(obj[key], "Array")) {
                    this._walkObjectToInclude(obj[key], dirpath);
                }
            }
        } else if (isInstanceOf(obj, "Array")) {
            for (let i = 0; i < obj.length; i++) {
                if (isInstanceOf(obj[i], "Object") || isInstanceOf(obj[i], "Array")) {
                    this._walkObjectToInclude(obj[i], dirpath);
                }
            }
        }
    }

    /**
     * Resolves `include` statements.
     * @param filePath {String} file to parse.
     */
    _parseIncludeStatements(filePath) {
        let data = safeParseJSON(fs.readFileSync(filePath, 'utf8'));
        if (OBJECT_ONLY && !isInstanceOf(data, "Object")) {
            throw "The JSON file being included should always be a dict rather than a list";
        }
        this._walkObjectToInclude(data, path.dirname(filePath));
        return data;
    }


    /**
     * Resolves includes and references
     * @param filePath {String} file to parse.
     */
    resolve(filePath) {
        const parsed = this._parseIncludeStatements(filePath);
        const dirName = path.dirname(filePath);
        let dereferenced = deref(parsed, {baseFolder: dirName});
        // handle circular references and use non-dereferenced source
        if ((dereferenced instanceof Error) && (dereferenced.message === "Circular self reference")) {
            dereferenced = parsed;
        }
        return dereferenced;
    }

    /**
     * Resolves includes and references for all the JSON files inside a given directory.
     * @param dirPath {String} directory to parse.
     */
    resolveDir(dirPath) {
        const data = [];
        file.walkSync(dirPath, (dirPath_, dirs_, files_) => {
            files_.forEach(file_ => data.push(this.resolve(path.join(dirPath_, file_))));
        });
        return data;
    }

}
