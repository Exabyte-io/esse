/* eslint-disable no-restricted-syntax */
/* eslint-disable no-await-in-loop */
import fs from "fs/promises";
import yaml from "js-yaml";
import path from "path";

const EnumType = new yaml.Type("!enum", {
    kind: "sequence",
    construct(data) {
        return {
            enum: data,
        };
    },
    represent(object) {
        return object.enum;
    },
});

const SCHEMA_DIR = "../../../schema/";
const yamlSchema = yaml.DEFAULT_SCHEMA.extend([EnumType]);

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
    if (path.extname(filePath) !== ".yml") {
        return;
    }

    const outFilename = `${path.basename(filePath, ".yml")}.json`;
    const dirname = path.dirname(filePath);
    const fileContents = await fs.readFile(filePath);
    const enumObj = yaml.load(fileContents, { schema: yamlSchema });
    await fs.writeFile(
        `${path.join(dirname, outFilename)}`,
        `${JSON.stringify(enumObj, null, 4)}\n`,
    );
});
