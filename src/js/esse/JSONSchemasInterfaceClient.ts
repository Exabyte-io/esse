import { SchemaObject } from "ajv";
import fs from "fs";
import path from "path";

import { walkDirSync } from "../scripts/utils";
import { JSONInterfaceQuery } from "./types";

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

export default class JSONSchemasInterfaceClient {
    static schemasCache = new Map<string, SchemaObject>();

    static addSchemaObject(schema: SchemaObject) {
        if (schema.$id) {
            this.schemasCache.set(schema.$id, schema);
        }
    }

    static schemaById(schemaId: string) {
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
    static matchSchema(query: JSONInterfaceQuery) {
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
}
