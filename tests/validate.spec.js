import chai from "chai";
import Ajv from "ajv";
import url from "url";
import path from "path";
import {walkSync} from "file";
import fs from "fs";
import _ from "lodash";
import prettyjson from "prettyjson";

const namespace = 'https://exabyte.io/schemas/';
const schemasDir = path.resolve(__dirname, '../schema');

const schemas = [];

const INCLUDE_KEY = '...';
const INCLUDE_VALUE_PATTERN = /^include\((.+)\)$/;
var JSON_INCLUDE_CACHE = {};

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
    console.log('Reading from', dirpath, filename);
    var filepath = path.join(dirpath, filename),
        json = fs.readFileSync(filepath, 'utf8'),
        d = JSON.parse(json);
    if (isInclude) {
        if (! isInstance(d, "Object")) {
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
 * Traverse schemas directory and build array of objects in following format:
 * {
 *     dir: 'primitive/',
 *     filename: '2d_data.json',
 *     content: {<schema definition is here>},
 *     example: {<example for the schema, optional>}
 * }
 * see walkSync docs for more info.
 */
walkSync(schemasDir, function (dirPath, dirs, files) {
    files.forEach(f=> {
        const filePath = path.join(dirPath, f);
        const dir = _.trimStart(dirPath.replace(schemasDir, ''), '/');

        const schema = {
            filename: f,
            content: JSON.parse(fs.readFileSync(filePath, 'utf8')),
            dir: dir === '' ? dir : dir + '/'
        };

        const exampleFile = filePath.replace('schema', 'example');
        if (fs.existsSync(exampleFile)) {
            JSON_INCLUDE_CACHE = {};
            var json = parseJSONInclude(path.dirname(exampleFile), path.basename(exampleFile));
            schema.example = json;
        }

        schemas.push(schema);
    });
});

// generate id for each schema and replace $refs values by appropriate id
schemas.forEach(schema=> {
    addIdProp(schema);
    replaceRefs(schema);
});

/**
 * Converts $ref values from "file:..." to "https://exabyte.io/schemas/..." (id).
 * @param schema
 */
function replaceRefs(schema) {
    var str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), url.resolve(namespace, schema.dir));
    schema.content = JSON.parse(str);
}

/**
 * Generates id from provided dir and filename. Id is equal to concatenation of: namespace, dir and filename.
 * @param dir {String}
 * @param filename {String}
 */
function generateId(dir, filename) {
    return url.resolve(namespace, path.join(dir, filename));
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

var ajv = Ajv({
    allErrors: true,
    schemas: schemas.map(schema=> schema.content)
});

describe('ajv', function () {
    it('works on primitive schema', function () {
        chai.expect(ajv.getSchema('https://exabyte.io/schemas/primitive/1d_data_series.json')).to.be.ok;
    });

    it('works on schema with external ref', function () {
        chai.expect(ajv.getSchema('https://exabyte.io/schemas/material/properties/structural/atomic_coordinates.json')).to.be.ok;
    });
});

describe('example', function () {
    schemas.forEach(schema => {
        if (schema.example) {
            const fullPath = path.join(schema.dir, schema.filename);
            it(`${fullPath} should be valid`, function () {
                const validator = ajv.getSchema(schema.id);
                const valid = validator(schema.example);
                if (!valid) {
                    console.log(JSON.stringify(schema.example))
                    console.log(prettyjson.render(validator.errors));
                }
                chai.expect(valid).to.be.ok;
            });
        }
    });
});
