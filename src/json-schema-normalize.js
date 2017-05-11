var url = require("url");
var path = require("path");
var walkSync = require("file").walkSync;
var fs = require("fs");
var _ = require("lodash");

var SCHEMAS_DIR = "../schema";
var NAMESPACE = "https://exabyte.io/schemas/";
var INCLUDE_KEY = '...';
var INCLUDE_VALUE_PATTERN = /^include\((.+)\)$/;
var JSON_INCLUDE_CACHE = {};

/**
 * @summary Is the object is instance of the type check.
 *
 * @param object the object
 * @param type the type
 * @return {boolean} the result
 */
function isInstance(object, type) {
    return (Object.prototype.toString.call(object).slice(8, -1) === type) ? true : false;
}

/**
 * Gets the name of the included file.
 * @param value String containing include statement
 */
function getIncludeName(value) {
    if ((isInstance(value, "String")) &&
        (value.search(INCLUDE_VALUE_PATTERN) != -1)) {
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
        json = fs.readFileSync(filepath, 'utf8'),
        d = JSON.parse(json);
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

/**
 * @summary Traverse schemas directory and build array of objects in following format:
 * ```json
 * {
 *     dir: 'primitive/', // relative folder location
 *     filename: '2d_data.json', // filename os the schema file
 *     content: {<schema definition is here>}, // schema content
 *     example: {<example for the schema, optional>}
 * }
 * ```
 * see walkSync docs for more info.
 *
 * @param withExamples if `true` add example object to the map
 */
function getNormalizedSchemas(withExamples) {

    var schemasDir = path.resolve(__dirname, SCHEMAS_DIR);
    var schemas = [];

    walkSync(schemasDir, function (dirPath, dirs, files) {
        files.forEach(function (f) {
            var filePath = path.join(dirPath, f);
            var dir = _.trimStart(dirPath.replace(schemasDir, ''), '/');

            var schema = {
                filename: f,
                content: JSON.parse(fs.readFileSync(filePath, 'utf8')),
                dir: dir === '' ? dir : dir + '/'
            };

            if (withExamples) {
                var exampleFile = filePath.replace('schema', 'example');
                if (fs.existsSync(exampleFile)) {
                    JSON_INCLUDE_CACHE = {};
                    schema.example = parseJSONInclude(path.dirname(exampleFile), path.basename(exampleFile));
                }
            }
            schemas.push(schema);
        });
    });

    // generate id for each schema and replace $refs values by appropriate id
    schemas.forEach(function (schema) {
        addIdProp(schema);
        replaceRefs(schema);
    });

    return schemas;
}

/**
 * Converts $ref values from "file:..." to "https://exabyte.io/schemas/..." (id).
 * @param schema
 */
function replaceRefs(schema) {
    var str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), url.resolve(NAMESPACE, schema.dir));
    schema.content = JSON.parse(str);
}

/**
 * Generates id from provided dir and filename. Id is equal to concatenation of: namespace, dir and filename.
 * @param dir {String}
 * @param filename {String}
 */
function generateId(dir, filename) {
    return url.resolve(NAMESPACE, path.join(dir, filename));
}

/**
 * Adds id property to schema if it does not exists.
 * @param schema
 */
function addIdProp(schema) {
    if (!schema.content.id) {
        schema.content.id = generateId(schema.dir, schema.filename);
    }
    schema.id = schema.content.id;
}

module.exports = {
    getNormalizedSchemas: getNormalizedSchemas,
    parseJSONInclude: parseJSONInclude
};

