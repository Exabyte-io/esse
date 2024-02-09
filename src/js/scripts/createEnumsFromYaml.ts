import fs from "fs";
import yaml from "js-yaml";
import lodash from "lodash";
import path from "path";

import { walkDir } from "../utils/filesystem";

/**
 *  We use YAML files to list enum options which need to be reused in other places.
 *  The format also allows for additional context to the often shortened slugs, such as DOIs.
 *  This script converts the enum lists to JSON files for $ref blocks inside schemas.
 */

const SCHEMA_DIR = "../../../schema/";

walkDir(SCHEMA_DIR, (filePath: string) => {
    if (path.extname(filePath) !== ".yml") {
        return;
    }
    console.log(filePath);
    const outFilename = `${path.basename(filePath, ".yml")}.json`;
    const dirname = path.dirname(filePath);

    const fileContents = fs.readFileSync(filePath);
    const obj = yaml.load(fileContents.toString()) as object;
    const enumObj = lodash.mapValues(obj, (value) => ({ enum: value }));

    fs.writeFileSync(`${path.join(dirname, outFilename)}`, `${JSON.stringify(enumObj, null, 4)}\n`);
});
