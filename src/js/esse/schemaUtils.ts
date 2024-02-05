// @ts-nocheck
import { SchemaObject } from "ajv";

export function mapObjectDeep(object: unknown, mapValue: (prop: unknown) => unknown): object {
    if (typeof object !== "object" || object === null) {
        return object;
    }

    if (Array.isArray(object)) {
        return object.map(
            (innerValue) => mapValue(innerValue) || mapObjectDeep(innerValue, mapValue),
        );
    }

    const entries = Object.entries(object).map(([key, value]) => {
        const res = mapValue(value);

        return [key, res === undefined ? mapObjectDeep(value, mapValue) : res];
    });

    return Object.fromEntries(entries);
}

export function walkSchema(object: unknown, callback: (prop: unknown) => unknown): object {
    if (Array.isArray(object)) {
        return object.map((item) => walkSchema(item, callback));
    }

    const cleanObject = callback(object);

    if (typeof cleanObject !== "object" || cleanObject === null) {
        return cleanObject;
    }

    const entries = Object.entries(cleanObject).map(([key, value]) => {
        return [key, walkSchema(value, callback)];
    });

    return Object.fromEntries(entries);
}

export function addAdditionalPropertiesToSchema(
    schema: SchemaObject,
    additionalProperties = false,
) {
    // @ts-ignore
    return walkSchema(schema, (object) => {
        if (typeof object !== "object") {
            return object;
        }

        const schema = object;

        if (schema.type === "object" && schema.properties && !("additionalProperties" in schema)) {
            return {
                ...schema,
                additionalProperties,
            };
        }

        return object;
    });
}

/**
 * By definition, the "compile" function will generate 2 schemas based on the following input:
 * {
 *    "title": "schema1",
 *    "type": "object",
 *    "properties": {
 *        "prop": {
 *             "title": "schema2",
 *             "type": "string"
 *        }
 *    }
 * }
 *
 * Result:
 * type schema2 = string;
 *
 * interface schema1 {
 *      prop: schema2
 * }
 *
 * To disable this behavior we need to remove "title" property from the "prop":
 *  {
 *    "title": "schema1",
 *    "type": "object",
 *    "properties": {
 *        "prop": {
 *             "type": "string"
 *        }
 *    }
 * }
 *
 * New result:
 * interface schema1 {
 *      prop: string;
 * }
 * @returns Clean schema
 */
export function cleanSchema<T = unknown>(object: T, clean = true): T {
    if (Array.isArray(object)) {
        return object.map((item) => cleanSchema(item)) as T;
    }

    if (typeof object !== "object" || object === null) {
        return object;
    }

    let cleanObject;

    if ((object as SchemaObject).title && clean) {
        // eslint-disable-next-line @typescript-eslint/no-unused-vars
        const { title, $schema, ...restObject } = object as SchemaObject;
        cleanObject = restObject;
    } else {
        cleanObject = object;
    }

    const entries = Object.entries(cleanObject).map(([key, value]) => {
        return [key, cleanSchema(value)];
    });

    return Object.fromEntries(entries);
}
