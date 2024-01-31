import fs from "fs/promises";
import path from "path";

import { patchSchema, walkDir } from "./utils.mjs";

const SCHEMA_DIR = "../../../schema/";

await walkDir(SCHEMA_DIR, async (filePath) => {
    if (path.extname(filePath) !== ".json") {
        return;
    }

    const fileContents = await fs.readFile(filePath);
    const schema = JSON.parse(fileContents);
    const schemaId = filePath.replace(SCHEMA_DIR, "").replace(".json", "").replace(/_/g, "-");
    let isFirstRun = true;
    const patchedSchema = patchSchema(schema, (object) => {
        if (isFirstRun) {
            isFirstRun = false;
            const { $id, ...restObject } = object;
            return {
                $id: schemaId,
                ...restObject,
                additionalProperties: object.additionalProperties || false,
            };
        }
        if (object.properties) {
            return { ...object, additionalProperties: object.additionalProperties || false };
        }
        return object;
    });

    const newContent = JSON.stringify(patchedSchema, null, 4);

    await fs.writeFile(filePath, `${newContent}\n`);
});
