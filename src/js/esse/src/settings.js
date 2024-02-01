import path from "path";
import { fileURLToPath } from "url";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

export const SCHEMAS_DIR = path.resolve(__dirname, "../../../../schema");
export const EXAMPLES_DIR = path.resolve(__dirname, "../../../../example");
export const PROPERTIES_MANIFEST_PATH = path.resolve(
    __dirname,
    "../../../../manifest/properties.yaml",
);
