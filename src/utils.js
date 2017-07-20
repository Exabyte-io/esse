import path from "path";
import fs from "fs";
import {LIB_DIR} from "./settings"

const INCLUDE_KEY = '...';
const INCLUDE_VALUE_PATTERN = /^include\((.+)\)$/;

function _isInstance(object, type) {
    return (Object.prototype.toString.call(object).slice(8, -1) === type);
}

// make sure that text representing arrays is parsed into Arrays
export function safeParseJSON(string) {
    let obj = JSON.parse(string);
    if (string[0] === "[") {
        return Object.keys(obj).map(function (key) { return obj[key]; });
    }
    return obj;
}

export class JSONSchemaResolver {

    constructor() {
        this.JSON_INCLUDE_CACHE = {};
    }

    _getIncludeFileName(value) {
        if ((_isInstance(value, "String")) &&
            (value.search(INCLUDE_VALUE_PATTERN) !== -1)) {
            return value.match(INCLUDE_VALUE_PATTERN)[1]
        }
    }

    _walkObjectToInclude(obj, dirpath) {
        if (_isInstance(obj, "Object")) {
            let isIncludeExp = false;
            if (INCLUDE_KEY in obj) {
                const includeName = this._getIncludeFileName(obj[INCLUDE_KEY])
                if (includeName) {
                    isIncludeExp = true;
                    delete obj[INCLUDE_KEY];
                    if (!(includeName in this.JSON_INCLUDE_CACHE)) {
                        const _f = path.join(dirpath, includeName);
                        this.JSON_INCLUDE_CACHE[includeName] = this.parseIncludeStatements(
                            path.dirname(_f), path.basename(_f), true
                        );
                        for (let attr in this.JSON_INCLUDE_CACHE[includeName]) {
                            obj[attr] = this.JSON_INCLUDE_CACHE[includeName][attr];
                        }
                    }
                }
            }
            if (isIncludeExp) {
                return;
            }
            for (let key in obj) {
                if (_isInstance(obj[key], "Object") || _isInstance(obj[key], "Array")) {
                    this._walkObjectToInclude(obj[key], dirpath);
                }
            }
        } else if (_isInstance(obj, "Array")) {
            for (let i = 0; i < obj.length; i++) {
                if (_isInstance(obj[i], "Object") || _isInstance(obj[i], "Array")) {
                    this._walkObjectToInclude(obj[i], dirpath);
                }
            }
        }
    }

    parseIncludeStatements(dirpath, filename, objectOnly, schemasList = []) {
        const relative = dirpath.replace(LIB_DIR, '').replace(/example|schema|/g, '').replace(/^\/+/, '');
        const filepath = path.join(dirpath, filename);
        let d;
        // either use passed list of schemas or read from disk
        if (schemasList.length) {
            const schema = schemasList.find((e) => {return e && e.dirpath === relative && e.filename === filename});
            d = schema.content;
        } else {
            const str = fs.readFileSync(filepath, 'utf8');
            d = safeParseJSON(str);
        }
        if (objectOnly) {
            if (!_isInstance(d, "Object")) {
                throw "The JSON file being included should always be a dict rather than a list";
            }
        }
        this._walkObjectToInclude(d, dirpath);
        return d;
    }
}

