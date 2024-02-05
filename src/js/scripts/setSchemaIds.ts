import fs from "fs";
import path from "path";

import { walkDirSync } from "./utils";

const SCHEMA_DIR = "../../../schema/";

walkDirSync(SCHEMA_DIR, (filePath) => {
    if (path.extname(filePath) !== ".json") {
        return;
    }

    const fileContents = fs.readFileSync(filePath);
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    const { $id, ...schema } = JSON.parse(fileContents.toString());
    const schemaId = filePath.replace(SCHEMA_DIR, "").replace(".json", "").replace(/_/g, "-");
    const newContent = JSON.stringify({ $id: schemaId, ...schema }, null, 4);

    fs.writeFileSync(filePath, `${newContent}\n`);
});
