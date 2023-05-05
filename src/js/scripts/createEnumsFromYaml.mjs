/* eslint-disable no-restricted-syntax */
import fs from "fs";
import yaml from "js-yaml";
import path from "path";
import lodash from "lodash";

/**
 *  We use YAML files to list enum options which need to be reused in other places.
 *  The format also allows for additional context to the often shortened slugs, such as DOIs.
 *  This script converts the enum lists to JSON files for $ref blocks inside schemas.
 */

const SCHEMA_DIR = "../../../schema/";

function walkDir(dir, callback) {
    const subDirs = fs.readdirSync(dir);

    for (const subDir of subDirs) {
        const itemPath = path.join(dir, subDir);
        const stat = fs.statSync(itemPath);

        if (stat.isDirectory()) {
            walkDir(itemPath, callback);
        } else {
            callback(itemPath);
        }
    }
}

walkDir(SCHEMA_DIR, (filePath) => {
    if (path.extname(filePath) !== ".yml") {
        return;
    }
    console.log(filePath);
    const outFilename = `${path.basename(filePath, ".yml")}.json`;
    const dirname = path.dirname(filePath);

    const fileContents = fs.readFileSync(filePath);
    const obj = yaml.load(fileContents);
    const enumObj = lodash.mapValues(obj, (value) => ({ enum: value }));

    fs.writeFileSync(
        `${path.join(dirname, outFilename)}`,
        `${JSON.stringify(enumObj, null, 4)}\n`,
    );
});
