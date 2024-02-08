import Ajv, { SchemaObject } from "ajv";

import { mapObjectDeep } from "../esse/schemaUtils";
import { JSONSchema } from "../esse/utils";

function addAdditionalPropertiesToSchema(schema: JSONSchema, additionalProperties = false) {
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
                unevaluatedProperties: false,
            };
        }
    });
}

const ajv = new Ajv({
    removeAdditional: true,
    strict: false, // TODO: adjust schemas and enable strict mode
    useDefaults: true,
    /**
     * discriminator fixes default values in oneOf
     * @see https://ajv.js.org/guide/modifying-data.html#assigning-defaults
     */
    discriminator: true,
});

export interface AnyObject {
    [key: string]: unknown;
}

export function getValidator(jsonSchema: SchemaObject) {
    const schemaKey = jsonSchema.$id as string;

    let validate = ajv.getSchema(schemaKey);

    if (!validate) {
        const patchedSchema = addAdditionalPropertiesToSchema(jsonSchema);
        ajv.addSchema(patchedSchema, schemaKey);
        validate = ajv.getSchema(schemaKey);
    }

    if (!validate) {
        throw new Error("JSONSchemasInterface AJV validator error");
    }

    return validate;
}

/**
 * Validates a given example against the schema.
 * @param example example to validate.
 * @param schema schema to validate the example with.
 * @returns whether example is valid.
 */
export function validate(data: AnyObject, jsonSchema: SchemaObject) {
    const validator = getValidator(jsonSchema);
    const isValid = validator(data);

    return {
        isValid,
        errors: validator.errors,
    };
}
