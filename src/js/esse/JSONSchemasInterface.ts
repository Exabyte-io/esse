import Ajv, { SchemaObject } from "ajv";
import fs from "fs";
import path from "path";

import { walkDirSync } from "../scripts/utils";
import { addAdditionalPropertiesToSchema } from "./schemaUtils";

type Query = { [key in keyof SchemaObject]: { $regex: string } };

export interface AnyObject {
    [key: string]: unknown;
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

export function readSchemaFolderSync(folderPath: string) {
    const schemas: SchemaObject[] = [];

    walkDirSync(folderPath, (filePath: string) => {
        if (path.extname(filePath) !== ".json") {
            return;
        }
        const schema = JSON.parse(fs.readFileSync(filePath).toString());
        schemas.push(schema);
    });

    return schemas;
}

export class JSONSchemasInterface {
    static schemaFolder = "./lib/js/schema";

    static schemasCache = new Map<string, SchemaObject>();

    static setSchemaFolder(schemaFolder: string) {
        if (this.schemaFolder !== schemaFolder) {
            this.schemaFolder = schemaFolder;
            this.readSchemaFolder();
        }
    }

    static readSchemaFolder() {
        const schemas = readSchemaFolderSync(this.schemaFolder);

        schemas.forEach((schema) => {
            if (schema.$id) {
                this.schemasCache.set(schema.$id, schema);
            }
        });
    }

    static schemaById(schemaId: string) {
        if (this.schemasCache.size === 0) {
            this.readSchemaFolder();
        }

        return this.schemasCache.get(schemaId);
    }

    /**
     * @example <caption>Search by $id regex</caption>
     * JSONSchemasInterface.matchSchema({
     *   $id: {
     *     $regex: 'software-application'
     *   }
     * })
     *
     * @example <caption>Search by $id and title regex</caption>
     * JSONSchemasInterface.matchSchema({
     *   $id: {
     *     $regex: 'software-application'
     *   },
     *   title: {
     *     $regex: 'application'
     *   }
     * })
     */
    static matchSchema(query: Query) {
        const searchFields = Object.keys(query) as Array<keyof typeof query>;

        return Array.from(this.schemasCache.values()).find((schema) => {
            return searchFields.every((field) => {
                const queryField = query[field];
                const schemaField = schema[field];

                if (!queryField || typeof schemaField !== "string") {
                    return;
                }

                return new RegExp(queryField.$regex).test(schemaField);
            });
        });
    }

    private static getAjvValidator(jsonSchema: SchemaObject) {
        const schemaKey = jsonSchema.$id as string;

        let validate = ajv.getSchema(schemaKey);

        if (!validate) {
            ajv.addSchema(addAdditionalPropertiesToSchema(jsonSchema), schemaKey);
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
    static validate(data: AnyObject, jsonSchema: SchemaObject) {
        const validator = this.getAjvValidator(jsonSchema);
        const isValid = validator(data);

        return {
            isValid,
            errors: validator.errors,
        };
    }
}
