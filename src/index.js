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

function writeFileToLibDir(content, relativePath) {
    const newPath = path.join(LIB_DIR, relativePath);
    if (!fs.existsSync(path.dirname(newPath))) {
        fs.mkdirSync(path.dirname(newPath));
    }
    fs.writeFileSync(newPath, JSON.stringify(content));
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

function getRawSchemasWithIncludeStatements() {
    const schemas = [];
    file.walkSync(SCHEMAS_DIR, function (dirPath, dirs, files) {
        files.forEach(function (f) {
            let schema = readSchemaFromPath(path.join(dirPath, f));
            schemas.push(replaceFileMentions(schema));
        });
    });
    return schemas;
}

const COMPILED_SCHEMAS = [],
    RAW_SCHEMAS = []; // raw schemas still contain `$ref` tags

const _schemas = getRawSchemasWithIncludeStatements();

_schemas.forEach((s, i, l) => {
    l[i] = parseIncludeStatements(s.dirpath, s.filename, false, l);
});

_schemas.forEach((s, i, l) => {
    const dirname = path.dirname(path.join(LIB_DIR, s.content.id));
    RAW_SCHEMAS.push(s.content);
    COMPILED_SCHEMAS.push(
        deref(s.content, {
            baseFolder: dirname
        })
    );
});

if (DEBUG) console.log(JSON.stringify(COMPILED_SCHEMAS, null, '\t'));

export const schemas = COMPILED_SCHEMAS;
export const rawSchemas = RAW_SCHEMAS;

export function getSchemaById(id, useRaw=False) {
    const s = raw ? RAW_SCHEMAS : COMPILED_SCHEMAS;
    return s.find(function (schema) {
        return schema.id === id;
    })
}

export function getSchemaByIdBasename(subId, useRaw=false) {
    const s = useRaw ? RAW_SCHEMAS : COMPILED_SCHEMAS;
    return s.find(function (schema) {
        return schema.id.split('/').reverse()[0] == subId;
    })
}
