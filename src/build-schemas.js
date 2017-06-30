var fs = require('file-system');
var path = require('path');
var url = require('url');
var _ = require("lodash");

var schemasDir = __dirname + '/../schema'; // relative path of schemas directory
var compiledSchemasPath = __dirname + '/schemas-list.json'; // relative path of compiled schemas

var allFiles = [];
var schemasList = [];

function traverseFileSystem(currentPath) {
    const files = fs.readdirSync(currentPath);
    for (var i in files) {
        const currentFile = currentPath + '/' + files[i];
        const stats = fs.statSync(currentFile);
        if (stats.isFile()) {
            allFiles.push(currentFile);
        }
        else if (stats.isDirectory()) {
            traverseFileSystem(currentFile);
        }
    }
}

var NAMESPACE = "https://exabyte.io/schemas/";

function replaceReferences(schema) {
    var str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), 'https:/schemas/' + schema.dir + '/');
    schema.content = JSON.parse(str);
}

function addIdProperty(schema) {
    if (!schema.content.id) {
        schema.content.id = url.resolve(NAMESPACE, path.join(schema.dir, schema.filename));
    }
    schema.id = schema.content.id;
}

function dereference(schema, schemasList) {
    var str = JSON.stringify(schema);
    if (str.indexOf('$ref') > 0) {
        // substitute globally ("$ref": "http:/schema/<path>") with resolved schema
        var regex = /[\"|\']\$ref[\'|\"]\s*\:\s*[\"|\']https\:\/schemas\/([^\"\'\s]+)[\'|\"]/g;
        str = str.replace(regex, function (match, subSchemaPath) {
            var s = _.find(schemasList, function (ss) {
                var subPath = path.resolve(subSchemaPath);
                var ssPath = ss.dir + '/' + ss.filename;
                console.log(ss.dir, ss.filename);
                return subPath.endsWith(ssPath);
            });
            if (!s) {
                throw new Error('Cannot find schema in schemas list: ' + match);
            }
            // replace de-referenced json, removing first '{' and last '}' character
            return JSON.stringify(dereference(s.content, schemasList)).replace(/^\s*\{\s*/, '').replace(/\s*\}\s*$/, '');
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

    return dereference(schema.content, schemasList);
}


traverseFileSystem(schemasDir);

// read all schemas
allFiles.forEach((el) => {
    schemasList.push(
        {
            dir: path.dirname(el.replace(schemasDir, '')),
            content: require(el),
            filename: path.basename(el)
        });
});

// prepare and dereference schemas
schemasList.forEach((schema, index, list) => {
    addIdProperty(schema);
    replaceReferences(schema);
    list[index] = dereference(schema, list);
});

fs.writeFileSync(compiledSchemasPath, JSON.stringify(schemasList));


