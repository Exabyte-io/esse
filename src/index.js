import fs from "fs";
import url from "url";
import path from "path";
import file from "file";
import lodash from "lodash";
import deref from "json-schema-deref-sync";
import {NAMESPACE, LIB_DIR, SCHEMAS_DIR, OMIT_SCHEMA_KEY, DEBUG} from "./settings";
import {parseIncludeStatements} from "./utils";

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
    const _dir = lodash.trimStart(path.dirname(p).replace(SCHEMAS_DIR, ''), '/'),
        schema = {
            filename: path.basename(p),
            content: JSON.parse(fs.readFileSync(p, 'utf8')),
            dirpath: _dir === '' ? _dir : _dir + '/',
        };
    // using relative path in filesystem wrt to the schemas root as id
    schema.content.id = path.join(schema.dirpath, schema.filename);
    if (OMIT_SCHEMA_KEY) {
        schema.content = lodash.omit(schema.content, "$schema");
    }
    return schema;
}

function getRawSchemas() {
    const schemas = [];
    file.walkSync(SCHEMAS_DIR, function (dirPath, dirs, files) {
        files.forEach(function (f) {
            let schema = readSchemaFromPath(path.join(dirPath, f));
            schemas.push(replaceFileMentions(schema));
        });
    });
    return schemas;
}

const ALL_SCHEMAS = [];
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

export const schemas = ALL_SCHEMAS;

export function getSchemaById(path) {
    return ALL_SCHEMAS.find(function (schema) {
        return schema.id === path;
    })
}
