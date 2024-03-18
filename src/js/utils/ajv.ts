import type { FuncKeywordDefinition } from "ajv";
import Ajv, { SchemaObject } from "ajv";
import { AnyValidateFunction } from "ajv/dist/core";

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

function addTransformDateKeywordToSchema(schema: JSONSchema) {
    return mapObjectDeep(schema, (object) => {
        const localSchema = object as JSONSchema;

        if (
            typeof object === "object" &&
            localSchema?.type === "string" &&
            localSchema?.format === "date-time"
        ) {
            // eslint-disable-next-line @typescript-eslint/no-unused-vars
            const { type, format, ...restSchema } = localSchema;
            return {
                ...restSchema,
                transformDate: true,
            };
        }
    });
}

/**
 * This function defines custom transformDate AJV keyword
 *
 * Problem:
 * There's no way to define dates except as a string with "date-time" or "date" format in JsonSchema:
 * {
 *    type: "string",
 *    format: "date-time"
 * }
 * But we need the Date object to be stored in databases or to perform correct date comparison
 *
 * Solution:
 * The keyword will convert string with "date-time" format to Date object
 */
function transformDateKeyword(): FuncKeywordDefinition {
    return {
        keyword: "transformDate",
        schemaType: "boolean",
        modifying: true,
        validate(transformDate, date, _metadata, dataCxt) {
            if (transformDate && dataCxt && typeof date === "string") {
                const dateObject = new Date(date);
                if (dateObject.toString() === "Invalid Date") {
                    return false;
                }

                dataCxt.parentData[dataCxt.parentDataProperty] = dateObject;
            }

            return true;
        },
    };
}

const ajvConfig = {
    strict: false, // TODO: adjust schemas and enable strict mode
    useDefaults: true,
    /**
     * discriminator fixes default values in oneOf
     * @see https://ajv.js.org/guide/modifying-data.html#assigning-defaults
     */
    discriminator: true,
    keywords: [transformDateKeyword()],
};

const ajvValidator = new Ajv({ ...ajvConfig });
const ajvValidatorAndCleaner = new Ajv({ ...ajvConfig, removeAdditional: true });
const ajvValidatorAndCleanerWithCoercingTypes = new Ajv({
    ...ajvConfig,
    removeAdditional: true,
    coerceTypes: true,
});

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
        // replace "date-time" format with "transformDate" keyword
        const patchedSchema = addTransformDateKeywordToSchema(
            // properties that were not defined in schema will be ignored when clean = false
            clean ? addAdditionalPropertiesToSchema(jsonSchema) : jsonSchema,
        );
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
