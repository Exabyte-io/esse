import fs from "fs";
import { compile } from "json-schema-to-typescript";

import { walkDir } from "./utils.mjs";

/**
 * By definition, the "compile" function will generate 2 schemas based on the following input:
 * {
 *    "title": "schema1",
 *    "type": "object",
 *    "properties": {
 *        "prop": {
 *             "title": "schema2",
 *             "type": "string"
 *        }
 *    }
 * }
 *
 * Result:
 * type schema2 = string;
 *
 * interface schema1 {
 *      prop: schema2
 * }
 *
 * To disable this behavior we need to remove "title" property from the "prop":
 *  {
 *    "title": "schema1",
 *    "type": "object",
 *    "properties": {
 *        "prop": {
 *             "type": "string"
 *        }
 *    }
 * }
 *
 * New result:
 * interface schema1 {
 *      prop: string;
 * }
 * @param {any} object
 * @param {Boolean} clean
 * @returns Clean schema
 */
function cleanSchema(object, clean = true) {
    if (Array.isArray(object)) {
        return object.map((item) => cleanSchema(item));
    }

    if (typeof object !== "object") {
        return object;
    }

    let cleanObject;

    if (object.title && clean) {
        const { title, $schema, ...restObject } = object;
        cleanObject = restObject;
    } else {
        cleanObject = object;
    }

    const entries = Object.entries(cleanObject).map(([key, value]) => {
        return [key, cleanSchema(value)];
    });

    return Object.fromEntries(entries);
}

export default async function compileTS(schemaPath, savePath) {
    try {
        await fs.promises.unlink(savePath);
    } catch (err) {
        console.log("File with types not exists");
    }

    await walkDir(schemaPath, async (filePath) => {
        const data = await fs.promises.readFile(filePath, "utf8");
        const schema = cleanSchema(JSON.parse(data), false);

        console.log(filePath);

        const compiledSchema = await compile(schema, schema.title, {
            unreachableDefinitions: true,
            additionalProperties: false,
            bannerComment: `/** Schema ${filePath} */`,
        });

        await fs.promises.appendFile(savePath, `${compiledSchema} \n`, { flag: "a+" });
    });
}
