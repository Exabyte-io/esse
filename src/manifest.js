import fs from "fs-extra";
import path from "path";
import {MANIFEST_DIR} from "./settings";

const MANIFEST_LIST = ["models.json", "properties.json"];

export default MANIFEST_LIST.map(file=> {
    const filePath = path.join(MANIFEST_DIR, file);
    return Object.assign(JSON.parse(fs.readFileSync(filePath)), {id: file, path: filePath});
});
