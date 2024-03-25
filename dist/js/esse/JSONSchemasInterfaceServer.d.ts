import { SchemaObject } from "ajv";
import JSONSchemasInterface from "./JSONSchemasInterface";
export declare function readSchemaFolderSync(folderPath: string): SchemaObject[];
export default class JSONSchemasInterfaceServer extends JSONSchemasInterface {
    static schemaFolder: string;
    static setSchemaFolder(schemaFolder: string): void;
    static readSchemaFolder(): void;
    static getSchemaById(schemaId: string): import("json-schema").JSONSchema7 | undefined;
}
