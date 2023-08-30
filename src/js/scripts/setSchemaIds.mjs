/* eslint-disable no-restricted-syntax */
/* eslint-disable no-await-in-loop */
import fs from "fs/promises";
import path from "path";

const SCHEMA_DIR = "../../../schema/";

async function walkDir(dir, callback) {
    const subDirs = await fs.readdir(dir);

    for (const subDir of subDirs) {
        const itemPath = path.join(dir, subDir);
        const stat = await fs.stat(itemPath);

        if (stat.isDirectory()) {
            await walkDir(itemPath, callback);
        } else {
            await callback(itemPath);
        }
    }
}

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
