import type {} from "ajv"; // @see https://github.com/microsoft/TypeScript/issues/47663
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

const ajvConfig = {
    strict: false, // TODO: adjust schemas and enable strict mode
    useDefaults: true,
    /**
     * discriminator fixes default values in oneOf
     * @see https://ajv.js.org/guide/modifying-data.html#assigning-defaults
     */
    discriminator: true,
};

const ajvValidator = new Ajv({ ...ajvConfig });
const ajvValidatorAndCleaner = new Ajv({ ...ajvConfig, removeAdditional: true });

export interface AnyObject {
    [key: string]: unknown;
}

export function getValidator(jsonSchema: SchemaObject, clean = false) {
    const schemaKey = jsonSchema.$id as string;

    const ajv = clean ? ajvValidatorAndCleaner : ajvValidator;

    let validate = ajv.getSchema(schemaKey);

    if (!validate) {
        // properties that were not defined in schema will be ignored when clean = false
        const patchedSchema = clean ? addAdditionalPropertiesToSchema(jsonSchema) : jsonSchema;
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

/**
 * Validates and clean a given example against the schema
 * @param example example to validate.
 * @param schema schema to validate the example with.
 * @returns whether example is valid.
 */
export function validateAndClean(data: AnyObject, jsonSchema: SchemaObject) {
    const validator = getValidator(jsonSchema, true);
    const isValid = validator(data);

    return {
        isValid,
        errors: validator.errors,
    };
}
