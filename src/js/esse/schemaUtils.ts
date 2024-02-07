import { SchemaObject } from "ajv";

import { JSONSchema, JSONSchemaDefinition } from "./utils";

export type MapSchema = (prop: JSONSchemaDefinition) => JSONSchemaDefinition | undefined;

export function mapObjectDeep(object: JSONSchemaDefinition, mapValue: MapSchema): JSONSchema;

export function mapObjectDeep(object: JSONSchemaDefinition[], mapValue: MapSchema): JSONSchema[];

export function mapObjectDeep(
    object: JSONSchemaDefinition | JSONSchemaDefinition[],
    mapValue: MapSchema,
): JSONSchemaDefinition | JSONSchemaDefinition[] {
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

export function addAdditionalPropertiesToSchema(schema: JSONSchema, additionalProperties = false) {
    return mapObjectDeep(schema, (object) => {
        const schema = object as JSONSchema;

        if (
            typeof object === "object" &&
            schema?.type === "object" &&
            schema?.properties &&
            !("additionalProperties" in schema)
        ) {
            return {
                ...schema,
                additionalProperties,
            };
        }
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
export function cleanSchema(schema: JSONSchema) {
    let firstRun = true;

    return mapObjectDeep(schema, (object) => {
        if (typeof object === "object" && object?.title && firstRun) {
            firstRun = false;
            // eslint-disable-next-line @typescript-eslint/no-unused-vars
            const { title, $schema, ...restObject } = object as SchemaObject;
            return restObject;
        }
    });
}
