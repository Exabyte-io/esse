import path from "path";

const INCLUDE_KEY = '...';
const INCLUDE_VALUE_PATTERN = /^include\((.+)\)$/;
const JSON_INCLUDE_CACHE = {};

function isInstance(object, type) {
    return (Object.prototype.toString.call(object).slice(8, -1) === type);
}

function getIncludeFileName(value) {
    if ((isInstance(value, "String")) &&
        (value.search(INCLUDE_VALUE_PATTERN) !== -1)) {
        return value.match(INCLUDE_VALUE_PATTERN)[1]
    }
}

function _walkObjectToInclude(obj, dirpath) {
    if (isInstance(obj, "Object")) {
        var isIncludeExp = false;
        if (INCLUDE_KEY in obj) {
            var includeName = getIncludeFileName(obj[INCLUDE_KEY])
            if (includeName) {
                isIncludeExp = true
                delete obj[INCLUDE_KEY]
                if (!(includeName in JSON_INCLUDE_CACHE)) {
                    var _f = path.join(dirpath, includeName);
                    JSON_INCLUDE_CACHE[includeName] = parseIncludeStatements(path.dirname(_f), path.basename(_f), true);
                    for (var attr in JSON_INCLUDE_CACHE[includeName]) {
                        obj[attr] = JSON_INCLUDE_CACHE[includeName][attr];
                    }
                }
            }
        }
        if (isIncludeExp) {
            return;
        }
        for (var key in obj) {
            if (isInstance(obj[key], "Object") || isInstance(obj[key], "Array")) {
                _walkObjectToInclude(obj[key], dirpath);
            }
        }
    } else if (isInstance(obj, "Array")) {
        for (var i = 0; i < obj.length; i++) {
            if (isInstance(obj[i], "Object") || isInstance(obj[i], "Array")) {
                _walkObjectToInclude(obj[i], dirpath);
            }
        }
    }
}

export function parseIncludeStatements(dirpath, filename, objectOnly, schemasList) {
    const filepath = path.join(dirpath, filename),
        // either use passed list of schemas or read from disk
        d = schemasList.length ? schemasList.find((el) => {
            return el.dirpath === dirpath && el.filename === filename
        }) : fs.readFileSync(filepath, 'utf8');
    if (objectOnly) {
        if (!isInstance(d, "Object")) {
            throw "The JSON file being included should always be a dict rather than a list";
        }
    }
    _walkObjectToInclude(d, dirpath);
    return d;
}
