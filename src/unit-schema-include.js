var url = require("url");
var path = require("path");
var _ = require("lodash");

var schemasList = [
    {
        dir: 'job/model/method/workflow',
        content: require('../schema/job/model/method/workflow/unit.json'),
        filename: 'unit.json'
    },
    {
        dir: 'job/model/method/workflow/unit',
        content: require('../schema/job/model/method/workflow/unit/map.json'),
        filename: 'map.json'
    },
    {
        dir: 'job/model/method/workflow/unit',
        content: require('../schema/job/model/method/workflow/unit/execution.json'),
        filename: 'execution.json'
    },
    {
        dir: 'job/model/method/workflow/unit',
        content: require('../schema/job/model/method/workflow/unit/assignment.json'),
        filename: 'assignment.json'
    },
    {
        dir: 'job/model/method/workflow/unit',
        content: require('../schema/job/model/method/workflow/unit/condition.json'),
        filename: 'condition.json'
    },
    {
        dir: 'job/model/method/workflow/unit',
        content: require('../schema/job/model/method/workflow/unit/convergence.json'),
        filename: 'convergence.json'
    },
    {
        dir: 'job/model/method/workflow/unit',
        content: require('../schema/job/model/method/workflow/unit/_result.json'),
        filename: '_result.json'
    },
    {
        dir: 'job/model/method/workflow/unit',
        content: require('../schema/job/model/method/workflow/unit/_input.json'),
        filename: '_input.json'
    },
    {
        dir: 'job',
        content: require('../schema/job/app.json'),
        filename: 'app.json'
    },
    {
        dir: 'job',
        content: require('../schema/job/compute.json'),
        filename: 'compute.json'
    },
    {
        dir: 'job/software',
        content: require('../schema/job/software/espresso.json'),
        filename: 'espresso.json'
    },
    {
        dir: 'job/software',
        content: require('../schema/job/software/vasp.json'),
        filename: 'vasp.json'
    },
    {
        dir: 'job/software',
        content: require('../schema/job/software/arguments.json'),
        filename: 'arguments.json'
    },
    {
        dir: 'job/software/espresso',
        content: require('../schema/job/software/espresso/arguments.json'),
        filename: 'arguments.json'
    }
];

var templateList = [
    {
        type: 'execution-espresso',
        template: require('../example/job/model/method/workflow/unit/templates/execution-espresso.json')
    },
    {
        type: 'execution-vasp',
        template: require('../example/job/model/method/workflow/unit/templates/execution-vasp.json')
    },
    {
        type: 'assignment',
        template: require('../example/job/model/method/workflow/unit/templates/assignment.json')
    },
    {
        type: 'condition',
        template: require('../example/job/model/method/workflow/unit/templates/condition.json')
    }
];

var NAMESPACE = "https://exabyte.io/schemas/";

schemasList.forEach(schema => {
    addIdProp(schema);
    replaceRefs(schema);
});

/**
 * Converts $ref values from "file:..." to "exabyte:<relative path to referenced schema>/...".
 * @param schema
 */
function replaceRefs(schema) {
    var str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), 'https:/schemas/' + schema.dir + '/');
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

function deref(schema, schemasList) {
    var str = JSON.stringify(schema);
    if (str.indexOf('$ref') > 0) {
        // substitute globally ("$ref": "http:/schema/<path>") with resolved schema
        var regex = /[\"|\']\$ref[\'|\"]\s*\:\s*[\"|\']https\:\/schemas\/([^\"\'\s]+)[\'|\"]/g;
        str = str.replace(regex, function (match, subSchemaPath) {
            var s = _.find(schemasList, function (ss) {
                var subPath = path.resolve(subSchemaPath);
                var ssPath = ss.dir + '/' + ss.filename;
                return subPath.endsWith(ssPath);
            });
            if (!s) {
                throw new Error('Cannot find schema in schemas list: ' + match);
            }
            // replace de-referenced json, removing first '{' and last '}' character
            return JSON.stringify(deref(s.content, schemasList)).replace(/^\s*\{\s*/, '').replace(/\s*\}\s*$/, '');
        });
        return JSON.parse(str);
    }
    return schema;
}

function jsonInclude(schemaPath) {
    const schema = _.find(schemasList, function (schema) {
        return schema.dir + '/' + schema.filename === schemaPath;
    });
    if (!schema) {
        throw new Error('Cannot find schema for path: ' + schemaPath);
    }

    return deref(schema.content, schemasList);
}

module.exports = {
    jsonInclude: jsonInclude,
    schemasList: schemasList,
    templateList: templateList
};
