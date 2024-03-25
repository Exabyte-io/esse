import Ajv, { SchemaObject } from "ajv";
import { AnyValidateFunction } from "ajv/dist/core";
import addFormats from "ajv-formats";

import { mapObjectDeep } from "../esse/schemaUtils";
import { AnyObject } from "../esse/types";
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
const ajvValidatorAndCleanerWithCoercingTypes = new Ajv({
    ...ajvConfig,
    removeAdditional: true,
    coerceTypes: true,
});

addFormats(ajvValidator);
addFormats(ajvValidatorAndCleaner);
addFormats(ajvValidatorAndCleanerWithCoercingTypes);

interface AjvInstanceOptions {
    clean: boolean;
    coerceTypes: boolean;
}

function getAjvInstance({ clean, coerceTypes }: AjvInstanceOptions) {
    if (clean && coerceTypes) {
        return ajvValidatorAndCleanerWithCoercingTypes;
    }

    if (clean) {
        return ajvValidatorAndCleaner;
    }

    return ajvValidator;
}

export function getValidator(
    jsonSchema: SchemaObject,
    { clean, coerceTypes }: AjvInstanceOptions,
): AnyValidateFunction {
    const schemaKey = jsonSchema.$id as string;
    const ajv = getAjvInstance({ clean, coerceTypes });
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
    const validator = getValidator(jsonSchema, { clean: false, coerceTypes: false });
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
export function validateAndClean(
    data: AnyObject,
    jsonSchema: SchemaObject,
    { coerceTypes = false },
) {
    const validator = getValidator(jsonSchema, { clean: true, coerceTypes });
    const isValid = validator(data);

    return {
        isValid,
        errors: validator.errors,
    };
}
