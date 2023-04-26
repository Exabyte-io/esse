/* eslint-disable no-restricted-syntax */
/* eslint-disable no-await-in-loop */
import fs from "fs/promises";
import yaml from "js-yaml";
import path from "path";
import lodash from "lodash";

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
    if (path.extname(filePath) !== ".yml" || path.extname(filePath) !== ".yaml") {
        return;
    }

    const outFilename = `${path.basename(filePath, ".yml")}.json`;
    const dirname = path.dirname(filePath);

    const fileContents = await fs.readFile(filePath);
    const obj = yaml.load(fileContents);
    const enumObj = lodash.mapValues(obj, (value) => ({ enum: value }));

    await fs.writeFile(
        `${path.join(dirname, outFilename)}`,
        `${JSON.stringify(enumObj, null, 4)}\n`,
    );
});
