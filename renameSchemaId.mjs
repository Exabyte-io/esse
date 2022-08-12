/* eslint-disable no-restricted-syntax */
/* eslint-disable no-await-in-loop */
import fs from "fs/promises";
import path from "path";

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

await walkDir("schema", async (filePath) => {
    if (path.extname(filePath) !== ".json") {
        return;
    }

    const fileContents = await fs.readFile(filePath);
    const { schemaId, ...schema } = JSON.parse(fileContents);

    if (!schemaId) {
        return;
    }

    const newSchemaId = schemaId.split("-").join("/");

    const newContent = JSON.stringify(
        {
            schemaId: newSchemaId,
            ...schema,
        },
        null,
        "    ",
    );

    await fs.writeFile(filePath, `${newContent}\n`);
});
