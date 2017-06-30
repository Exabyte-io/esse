var url = require("url");
var path = require("path");
var walkSync = require("file").walkSync;
var fs = require("fs");
var _ = require("lodash");
var deref = require('json-schema-deref-sync');
var Resolver = require('json-schema-ref-parser');

var SCHEMAS_DIR = "../schema";
var LIB_DIR = "../lib";
var NAMESPACE = "https://exabyte.io/schemas/";
var INCLUDE_KEY = '...';
var INCLUDE_VALUE_PATTERN = /^include\((.+)\)$/;
var JSON_INCLUDE_CACHE = {};
var ALL_SCHEMAS = [];

var schemasDir = path.resolve(__dirname, SCHEMAS_DIR);
var libDir = path.resolve(__dirname, LIB_DIR);

/**
 * @summary Is the object is instance of the type check.
 *
 * @param object the object
 * @param type the type
 * @return {boolean} the result
 */
function isInstance(object, type) {
    return (Object.prototype.toString.call(object).slice(8, -1) === type);
}

/**
 * Gets the name of the included file.
 * @param value String containing include statement
 */
function getIncludeName(value) {
    if ((isInstance(value, "String")) &&
        (value.search(INCLUDE_VALUE_PATTERN) !== -1)) {
        return value.match(INCLUDE_VALUE_PATTERN)[1]
    }
}

/**
 * Parses include statements inside a file.
 * @param dirpath Path where the include statements should be resolved in
 * @param filename Name of the json file to read in
 */
function parseJSONInclude(dirpath, filename, isInclude) {
    var filepath = path.join(dirpath, filename),
        // either use passed list of schemas or read from disk
        json = ALL_SCHEMAS.length ? ALL_SCHEMAS.find((el) => {
            return el.dirpath == dirpath && el.filename == filename
        }) : fs.readFileSync(filepath, 'utf8');
    d = json;
    if (isInclude) {
        if (!isInstance(d, "Object")) {
            throw "The JSON file being included should always be a dict rather than a list";
        }
    }
    walkThoughToInclude(d, dirpath);
    return d;
}

/**
 * Converts "...": "include(file.json)" values from to json content of "file.json".
 * @param obj JSON object containing include statements
 * @param dirpath Path where the include statements should be resolved in
 */
function walkThoughToInclude(obj, dirpath) {
    if (isInstance(obj, "Object")) {
        var isIncludeExp = false;
        if (INCLUDE_KEY in obj) {
            var includeName = getIncludeName(obj[INCLUDE_KEY])
            if (includeName) {
                isIncludeExp = true
                delete obj[INCLUDE_KEY]
                if (!(includeName in JSON_INCLUDE_CACHE)) {
                    var _f = path.join(dirpath, includeName);
                    JSON_INCLUDE_CACHE[includeName] = parseJSONInclude(path.dirname(_f), path.basename(_f), true);
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
            if (isInstance(obj[key], "Object") ||
                isInstance(obj[key], "Array")) {
                walkThoughToInclude(obj[key], dirpath);
            }
        }
    } else if (isInstance(obj, "Array")) {

        for (var i = 0; i < obj.length; i++) {
            if (isInstance(obj[i], "Object") ||
                isInstance(obj[i], "Array")) {
                walkThoughToInclude(obj[i], dirpath);
            }
        }
    }
}

function replaceFileRefsWithUrl(schema) {
    var str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), url.resolve(NAMESPACE, schema.dirpath));
    return Object.assign({}, schema, {content: JSON.parse(str)});
}

function replaceFileRefsWithNone(schema) {
    var str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), '');
    return Object.assign({}, schema, {content: JSON.parse(str)});
}

function addIdProperty(schema) {
    if (!schema.content.id) {
        schema.content.id = url.resolve(NAMESPACE, path.join(schema.dirpath, schema.filename));
    }
    return Object.assign({}, schema, {id: schema.content.id});
}

function getDereferencedSchemas(withExamples) {

    var schemas = [];

    walkSync(schemasDir, function (dirPath, dirs, files) {
        files.forEach(function (f) {
            var filePath = path.join(dirPath, f);
            var _dir = _.trimStart(dirPath.replace(schemasDir, ''), '/');

            var schema = {
                filename: f,
                content: JSON.parse(fs.readFileSync(filePath, 'utf8')),
                dirpath: _dir === '' ? _dir : _dir + '/'
            };

            schema = replaceFileRefsWithNone(schema)

            schemas.push(schema);

            var newDir = path.join(libDir, schema.dirpath);
            var newFil = path.join(newDir, schema.filename);
            if (!fs.existsSync(newDir)) {
                fs.mkdirSync(newDir);
            }
            fs.writeFileSync(newFil, JSON.stringify(schema));
        });
    });

    return schemas;
}

ALL_SCHEMAS = getDereferencedSchemas().slice(0, 1);

ALL_SCHEMAS.forEach((s, i, l) => {
    l[i] = parseJSONInclude(s.dirpath, s.filename);
});

ALL_SCHEMAS.forEach((s, i, l) => {
    l[i] = deref(s, {
        baseFolder: libDir
    });
//    Resolver.dereference(s, {
//        resolve: {
//            file: {
//                canRead: new RegExp(libDir, "i")
//            }
//        }
//    })
//        .then((schema) => {
//            l[i] = schema
//        })
//        .catch((err) => {
//            console.error(err)
//        })
});

console.log(JSON.stringify(ALL_SCHEMAS, null, '\t'));

module.exports = {
    getNormalizedSchemas: getDereferencedSchemas,
    parseJSONInclude: parseJSONInclude
};
