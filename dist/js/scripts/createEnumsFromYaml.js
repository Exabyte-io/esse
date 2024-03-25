"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
const fs_1 = __importDefault(require("fs"));
const js_yaml_1 = __importDefault(require("js-yaml"));
const lodash_1 = __importDefault(require("lodash"));
const path_1 = __importDefault(require("path"));
const filesystem_1 = require("../utils/filesystem");
/**
 *  We use YAML files to list enum options which need to be reused in other places.
 *  The format also allows for additional context to the often shortened slugs, such as DOIs.
 *  This script converts the enum lists to JSON files for $ref blocks inside schemas.
 */
const SCHEMA_DIR = "../../../schema/";
(0, filesystem_1.walkDir)(SCHEMA_DIR, (filePath) => {
    if (path_1.default.extname(filePath) !== ".yml") {
        return;
    }
    console.log(filePath);
    const outFilename = `${path_1.default.basename(filePath, ".yml")}.json`;
    const dirname = path_1.default.dirname(filePath);
    const fileContents = fs_1.default.readFileSync(filePath);
    const obj = js_yaml_1.default.load(fileContents.toString());
    const enumObj = lodash_1.default.mapValues(obj, (value) => ({ enum: value }));
    fs_1.default.writeFileSync(`${path_1.default.join(dirname, outFilename)}`, `${JSON.stringify(enumObj, null, 4)}\n`);
});
