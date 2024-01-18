/**
 * build_schemas uses node API to read all schemas from the FS
 * at build time and writes them out to a single schemas.js file for
 * downstream consumption to avoid FS calls in the browser.
 */
const fs = require("fs");
const path = require("path");
const mergeAllOf = require("json-schema-merge-allof");
const { ESSE } = require("./lib/js/esse");

const esse = new ESSE();
const { schemas, wrappedExamples } = esse;
const schema = esse.buildGlobalSchema();

fs.writeFileSync(
    "./schemas.js",
    "module.exports = {schemas: " + JSON.stringify(schemas) + "}",
    "utf8",
);

fs.writeFileSync("./schema.js", "module.exports = " + JSON.stringify(schema), "utf8");

if (process.env.BUILD_ASSETS !== "true") {
    process.exit(0);
}

const subfolder = process.env.BUILD_PATH || "./docs/js/";
schemas.forEach((s) => {
    if (process.env.SKIP_MERGE_ALLOF !== "true") {
        s = mergeAllOf(s, { resolvers: { defaultResolver: mergeAllOf.options.resolvers.title } });
    }
    const id_as_path = s.$id.replace("-", "_");
    const full_path = `${subfolder}/schema/${id_as_path}.json`;
    fs.mkdirSync(path.dirname(full_path), { recursive: true });
    fs.writeFileSync(full_path, JSON.stringify(s, null, 4), "utf8");
});
wrappedExamples.forEach((e) => {
    const id_as_path = e.path.replace("-", "_");
    const full_path = `${subfolder}/example/${id_as_path}.json`;
    fs.mkdirSync(path.dirname(full_path), { recursive: true });
    fs.writeFileSync(full_path, JSON.stringify(e.data, null, 4), "utf8");
});
