import fs from "fs";
import { compile, JSONSchema } from "json-schema-to-typescript";

import { cleanSchema } from "../esse/schemaUtils";
import { walkDir } from "./utils";

export default async function compileTS(schemaPath: string, savePath: string) {
    try {
        await fs.promises.unlink(savePath);
    } catch (err) {
        console.log("File with types not exists");
    }

    await walkDir(schemaPath, async (filePath) => {
        const data = await fs.promises.readFile(filePath, "utf8");
        const schema = cleanSchema(JSON.parse(data));

        console.log(filePath);

        const compiledSchema = await compile(schema as JSONSchema, schema.title || "", {
            unreachableDefinitions: true,
            additionalProperties: false,
            bannerComment: `/** Schema ${filePath} */`,
        });

        await fs.promises.appendFile(savePath, `${compiledSchema} \n`, { flag: "a+" });
    });
}
