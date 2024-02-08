import { SchemaObject } from "ajv";

export type JSONInterfaceQuery = { [key in keyof SchemaObject]: { $regex: string } };

export interface AnyObject {
    [key: string]: unknown;
}
