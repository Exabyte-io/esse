import { SchemaObject } from "ajv";
import { AnyValidateFunction } from "ajv/dist/core";
import { AnyObject } from "../esse/types";
interface AjvInstanceOptions {
    clean: boolean;
    coerceTypes: boolean;
}
export declare function getValidator(jsonSchema: SchemaObject, { clean, coerceTypes }: AjvInstanceOptions): AnyValidateFunction;
/**
 * Validates a given example against the schema.
 * @param example example to validate.
 * @param schema schema to validate the example with.
 * @returns whether example is valid.
 */
export declare function validate(data: AnyObject, jsonSchema: SchemaObject): {
    isValid: boolean | Promise<any>;
    errors: import("ajv").ErrorObject<string, Record<string, any>, unknown>[] | null | undefined;
};
/**
 * Validates and clean a given example against the schema
 * @param example example to validate.
 * @param schema schema to validate the example with.
 * @returns whether example is valid.
 */
export declare function validateAndClean(data: AnyObject, jsonSchema: SchemaObject, { coerceTypes }: {
    coerceTypes?: boolean | undefined;
}): {
    isValid: boolean | Promise<any>;
    errors: import("ajv").ErrorObject<string, Record<string, any>, unknown>[] | null | undefined;
};
export {};
