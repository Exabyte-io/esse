import {getRawJSONWithIncludeStatements,includeAndDereferenceJSONData} from "./deref";

const COMPILED_SCHEMAS = [],
    RAW_SCHEMAS = [], // raw schemas still contain `$ref` tags
    EXAMPLES = [];

const _examples = getRawJSONWithIncludeStatements(true);
const _schemas = getRawJSONWithIncludeStatements();

includeAndDereferenceJSONData({list: _examples, rawStore: EXAMPLES, example: true});
includeAndDereferenceJSONData({list: _schemas, compiledStore: COMPILED_SCHEMAS, rawStore: RAW_SCHEMAS});

// if (process.env.DEBUG) console.log(JSON.stringify(COMPILED_SCHEMAS, null, '\t'));

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
        return schema && schema.id && schema.id.split('/').reverse()[0] === subId;
    })
}

export function getExampleByIdBasename(subId) {
    return EXAMPLES.find(function (example) {
        return example && example.id && example.id.split('/').reverse()[0] === subId;
    })
}
