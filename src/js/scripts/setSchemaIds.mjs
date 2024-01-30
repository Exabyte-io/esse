import fs from "fs/promises";
import path from "path";

import { walkDir } from "./utils.mjs";

const SCHEMA_DIR = "../../../schema/";

await walkDir(SCHEMA_DIR, async (filePath) => {
    if (path.extname(filePath) !== ".json") {
        return;
    }

    const fileContents = await fs.readFile(filePath);
    const { $id: oldSchemaId, ...schema } = JSON.parse(fileContents);
    const schemaId = filePath.replace(SCHEMA_DIR, "").replace(".json", "").replace(/_/g, "-");
    const newContent = JSON.stringify({ $id: schemaId, ...schema }, null, 4);

    await fs.writeFile(filePath, `${newContent}\n`);
});
