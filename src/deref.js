import fs from "fs-extra";
import url from "url";
import path from "path";
import file from "file";
import lodash from "lodash";
import deref from "json-schema-deref-sync";
import {NAMESPACE, LIB_DIR, SCHEMAS_DIR, EXAMPLES_DIR, OMIT_SCHEMA_KEY, DEBUG} from "./settings";
import {JSONSchemaResolver, safeParseJSON} from "./utils";

function replaceFileMentions(schema, withUrl) {
    const suffix = withUrl ? url.resolve(NAMESPACE, schema.dirpath) : '',
        str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), suffix);
    return Object.assign({}, schema, {content: safeParseJSON(str)});
}

function writeFileToLibDir(content, relPath) {
    const fullPath = path.join(LIB_DIR, relPath);
    if (!fs.existsSync(path.dirname(fullPath))) {
        fs.ensureDirSync(path.dirname(fullPath));
    }
    fs.writeFileSync(fullPath, JSON.stringify(content));
}

function readJSONFromPath(p, example=false) {
    const d = example ? EXAMPLES_DIR : SCHEMAS_DIR,
        fileContent = fs.readFileSync(p, 'utf8'),
        data = {
            filename: path.basename(p),
            content: safeParseJSON(fileContent),
            dirpath: path.dirname(p).replace(`${d}`, '').replace(/^\//, ''),
        };
    // using relative path in filesystem wrt to the schema|example dir as id
    if (example) {
        // examples can be arrays -> add id to top level
        data.id = path.join(data.dirpath, data.filename);
    } else {
        // for schemas add `id` to `content` directly for future dereferencing
        data.content.id = path.join(data.dirpath, data.filename);
        if (OMIT_SCHEMA_KEY) data.content = lodash.omit(data.content, "$schema");
    }
    return data;
}

export function getRawJSONWithIncludeStatements(example=false) {
    const d = example ? EXAMPLES_DIR : SCHEMAS_DIR;
    const prefix = example ? "example" : "schema";
    const schemas = [];
    file.walkSync(d, function (dirPath, dirs, files) {
        files.forEach(function (f) {
            let data = readJSONFromPath(path.join(dirPath, f), example);
            data = replaceFileMentions(data);
            writeFileToLibDir(data.content, path.join(prefix, example ? data.id : data.content.id));
            schemas.push(data);
        });
    });
    return schemas;
}

export function includeAndDereferenceJSONData({list, compiledStore, rawStore=[], example=false}) {
    const prefix = example ? "example" : "schema";
    list.forEach((s, i, l) => {
        const dirname = path.dirname(path.join(LIB_DIR, prefix, example ? s.id : s.content.id));
        const IncludeResolver = new JSONSchemaResolver();
        const parsed = IncludeResolver.parseIncludeStatements(dirname, s.filename, false, rawStore, example);
        rawStore.push(example ? {content: parsed, id: s.id} : parsed);
    });

    if (compiledStore) {
        rawStore.forEach((el) => {
            const dirname = path.dirname(path.join(LIB_DIR, prefix, el.id));
            compiledStore.push(deref(el, {baseFolder: dirname}));
        });
    }
}
