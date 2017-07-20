import fs from "fs";
import fse from "fs-extra";
import url from "url";
import path from "path";
import file from "file";
import lodash from "lodash";
import deref from "json-schema-deref-sync";
import {NAMESPACE, LIB_DIR, SCHEMAS_DIR, EXAMPLES_DIR, OMIT_SCHEMA_KEY, DEBUG} from "./settings";
import {parseIncludeStatements} from "./utils";

function replaceFileMentions(schema, withUrl) {
    const suffix = withUrl ? url.resolve(NAMESPACE, schema.dirpath) : '',
        str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), suffix);
    return Object.assign({}, schema, {content: JSON.parse(str)});
}

function writeFileToLibDir(content, relPath) {
    const fullPath = path.join(LIB_DIR, relPath);
    if (!fs.existsSync(path.dirname(fullPath))) {
        fse.ensureDirSync(path.dirname(fullPath));
    }
    fs.writeFileSync(fullPath, JSON.stringify(content));
}

function readJSONFromPath(p, example=false) {
    const d = example ? EXAMPLES_DIR : SCHEMAS_DIR,
        schema = {
            filename: path.basename(p),
            content: JSON.parse(fs.readFileSync(p, 'utf8')),
            dirpath: path.dirname(p).replace(`${d}`, '').replace(/^\//, ''),
        };
    // using relative path in filesystem wrt to the schemas root as id
    schema.content.id = path.join(schema.dirpath, schema.filename);
    if (OMIT_SCHEMA_KEY) {
        schema.content = lodash.omit(schema.content, "$schema");
    }
    return schema;
}

function getRawJSONWithIncludeStatements(example=false) {
    const d = example ? EXAMPLES_DIR : SCHEMAS_DIR;
    const prefix = example ? "example" : "schema";
    const schemas = [];
    file.walkSync(d, function (dirPath, dirs, files) {
        files.forEach(function (f) {
            let schema = readJSONFromPath(path.join(dirPath, f), example);
            schema = replaceFileMentions(schema);
            writeFileToLibDir(schema.content, path.join(prefix, schema.content.id));
            schemas.push(schema);
        });
    });
    return schemas;
}

function dereferenceJSONData({list, compiledStore, rawStore, example=false}) {
    const prefix = example ? "example" : "schema";
    list.forEach((s, i, l) => {
        const dirname = path.dirname(path.join(LIB_DIR, prefix, s.content.id));
        const parsed = parseIncludeStatements(dirname, s.filename, false, l, example);
        console.log(i, parsed);
        l[i] = parsed;
    });

    list.forEach((s, i, l) => {
        console.log(s, i);
        const dirname = path.dirname(path.join(LIB_DIR, prefix, s.content.id));
        compiledStore.push(
            deref(s.content, {
                baseFolder: dirname
            })
        );
        if (rawStore) rawStore.push(s.content);
    });
}

const COMPILED_SCHEMAS = [],
    RAW_SCHEMAS = [], // raw schemas still contain `$ref` tags
    EXAMPLES = [];

const _schemas = getRawJSONWithIncludeStatements();
const _examples = getRawJSONWithIncludeStatements(true);

dereferenceJSONData({list: _examples, compiledStore: EXAMPLES, example: true});
dereferenceJSONData({list: _schemas, compiledStore: COMPILED_SCHEMAS, rawStore: RAW_SCHEMAS});

if (DEBUG) console.log(JSON.stringify(COMPILED_SCHEMAS, null, '\t'));

export const schemas = COMPILED_SCHEMAS;
export const rawSchemas = RAW_SCHEMAS;
export const examples = EXAMPLES;

export function getSchemaById(id, useRaw=False) {
    const s = raw ? RAW_SCHEMAS : COMPILED_SCHEMAS;
    return s.find(function (schema) {
        return schema.id === id;
    })
}

export function getSchemaByIdBasename(subId, useRaw=false) {
    const s = useRaw ? RAW_SCHEMAS : COMPILED_SCHEMAS;
    return s.find(function (schema) {
        return schema.id.split('/').reverse()[0] === subId;
    })
}

export function getExampleByIdBasename(subId) {
    return EXAMPLES.find(function (example) {
        return example.id.split('/').reverse()[0] === subId;
    })
}
