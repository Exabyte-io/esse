const url = require("url");
const path = require("path");
const walkSync = require("file").walkSync;
const fs = require("fs");
const _ = require("lodash");
const deref = require("json-schema-deref-sync");

/****************************************
 *              Set up                  *
 ****************************************/
const DEBUG = process.env['DEBUG'] || false;
const NAMESPACE = "https://exabyte.io/schemas/";

const INCLUDE_KEY = '...';
const INCLUDE_VALUE_PATTERN = /^include\((.+)\)$/;

const SCHEMAS_DIR = path.resolve(__dirname, "../schema");
const LIB_DIR = path.resolve(__dirname, "../lib");

const JSON_INCLUDE_CACHE = {};
const ALL_SCHEMAS = [];
const OMIT_SCHEMA_KEY = true;

/****************************************
 *     Schema manupulation tools        *
 ****************************************/
function isInstance(object, type) {
    return (Object.prototype.toString.call(object).slice(8, -1) === type);
}

function getIncludeFileName(value) {
    if ((isInstance(value, "String")) &&
        (value.search(INCLUDE_VALUE_PATTERN) !== -1)) {
        return value.match(INCLUDE_VALUE_PATTERN)[1]
    }
}

function parseIncludeStatements(dirpath, filename, objectOnly, schemasList) {
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

function replaceFileMentions(schema, withUrl) {
    const suffix = withUrl ? url.resolve(NAMESPACE, schema.dirpath) : '',
        str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), suffix);
    return Object.assign({}, schema, {content: JSON.parse(str)});
}

function writeFileToLibDir(schema) {
    const newDir = path.join(LIB_DIR, schema.dirpath),
        newFil = path.join(newDir, schema.filename);
    if (!fs.existsSync(newDir)) {
        fs.mkdirSync(newDir);
    }
    fs.writeFileSync(newFil, JSON.stringify(schema.content));
}

function readSchemaFromPath(p) {
    const _dir = _.trimStart(path.dirname(p).replace(SCHEMAS_DIR, ''), '/'),
        schema = {
            filename: path.basename(p),
            content: JSON.parse(fs.readFileSync(p, 'utf8')),
            dirpath: _dir === '' ? _dir : _dir + '/',
        };
    // using relative path in filesystem wrt to the schemas root as id
    schema.content.id = path.join(schema.dirpath, schema.filename);
    if (OMIT_SCHEMA_KEY) {
        schema.content = _.omit(schema.content, "$schema");
    }
    return schema;
}

function getRawSchemas() {
    const schemas = [];
    walkSync(SCHEMAS_DIR, function (dirPath, dirs, files) {
        files.forEach(function (f) {
            let schema = readSchemaFromPath(path.join(dirPath, f));
            schemas.push(replaceFileMentions(schema));
        });
    });
    return schemas;
}

/****************************************
 * Execute the logic and export schemas *
 ****************************************/

const _schemas = getRawSchemas();

_schemas.forEach((s, i, l) => {
    l[i] = parseIncludeStatements(s.dirpath, s.filename, false, l);
});

_schemas.forEach((s, i, l) => {
    const dirname = path.dirname(path.join(LIB_DIR, s.content.id));
    ALL_SCHEMAS.push(
        deref(s.content, {
            baseFolder: dirname
        })
    );
});

if (DEBUG) console.log(JSON.stringify(ALL_SCHEMAS, null, '\t'));

module.exports = {
    schemas: ALL_SCHEMAS,
    getSchemaByPath: function (path) {
        return ALL_SCHEMAS.find(function (schema) {
            return schema.path === path;
        })
    }
};
