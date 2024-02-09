import { SchemaObject } from "ajv";
import fs from "fs";
import path from "path";

import { walkDirSync } from "../utils/filesystem";
import JSONSchemasInterface from "./JSONSchemasInterface";

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

export default class JSONSchemasInterfaceServer extends JSONSchemasInterface {
    static schemaFolder = path.resolve(__dirname, "./../schema");

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

    static getSchemaById(schemaId: string) {
        if (this.schemasCache.size === 0) {
            this.readSchemaFolder();
        }

        return super.getSchemaById(schemaId);
    }
}
