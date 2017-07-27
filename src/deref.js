import fs from "fs-extra";
import url from "url";
import path from "path";
import file from "file";
import lodash from "lodash";
import deref from "json-schema-deref-sync";
import {NAMESPACE, LIB_DIR, SCHEMAS_DIR, EXAMPLES_DIR, OMIT_SCHEMA_KEY, DEBUG} from "./settings";
import {JSONSchemaResolver, safeParseJSON} from "./utils";

/**
 * Replaces mentions of `file:` inside schema references with empty string or url from namespace
 * @param {Object} schema object
 * @param {Boolean} withUrl whether to replace with URL
 */
function replaceFileMentions(schema, withUrl) {
    const suffix = withUrl ? url.resolve(NAMESPACE, schema.dirpath) : '',
        str = JSON.stringify(schema.content).replace(new RegExp('file:', 'g'), suffix);
    return Object.assign({}, schema, {content: safeParseJSON(str)});
}
/**
 * Writes textual `content` to LIB_DIR directory using relative path `relPath`
 * @param {String} content Textual content of the file
 * @param {String} relPath Relative path to the file within LIB_DIR
 */
function writeFileToLibDir(content, relPath) {
    const fullPath = path.join(LIB_DIR, relPath);
    if (!fs.existsSync(path.dirname(fullPath))) {
        fs.ensureDirSync(path.dirname(fullPath));
    }
    fs.writeFileSync(fullPath, JSON.stringify(content));
}

/**
 * Reads json file, either schema or example from path
 * @param {String} p Full path in filesystem
 * @param {Boolean} example Whether to read example file, otherwise assume schemas
 */
function readJSONFromPath(p, example = false) {
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

/**
 * Get JSON Data as its stored on disk, eg. raw without resolving references and include statements
 * @param {Boolean} example Whether to get examples, otherwise assume schemas
 */
export function getRawJSONWithIncludeStatements(example = false) {
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

/**
 * Reads json file, either schema or example from path
 * @param {Array} list List of raw schemas
 * @param {Array} compiledStore List used to store compiled data, modified as a result
 * @param {Array} rawStore List used to store raw data (without dereferencing, only resolving `include` statements)
 */
export function includeAndDereferenceJSONData({list, compiledStore, rawStore = [], example = false}) {
    const prefix = example ? "example" : "schema";
    list.forEach((s, i, l) => {
        const dirname = path.dirname(path.join(LIB_DIR, prefix, example ? s.id : s.content.id));
        const IncludeResolver = new JSONSchemaResolver();
        const parsed = IncludeResolver.parseIncludeStatements(dirname, s.filename, false, rawStore, example);
        rawStore.push(example ? {
            content: parsed,
            id: s.id
        } : parsed);
    });

    if (compiledStore) {
        rawStore.forEach((el) => {
            const dirname = path.dirname(path.join(LIB_DIR, prefix, el.id));
            let dereferenced = deref(el, {baseFolder: dirname});
            // handle circular references and use non-dereferenced source
            if ((dereferenced instanceof Error) && (dereferenced.message === "Circular self reference")) {
                dereferenced = el;
            }
            compiledStore.push(dereferenced);
        });
    }
}
