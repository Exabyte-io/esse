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

if (process.env.BUILD_DOCS !== "true") {
    process.exit(0);
}
schemas.forEach((schema) => {
    if (!process.env.SKIP_MERGE_ALLOF === "true") {
        schema = mergeAllOf(schema, {resolvers: {defaultResolver: mergeAllOf.options.resolvers.title}});
    }
    id_as_path = schema["$id"].replace("-", "_");
    full_path = `./docs/js/schema/${id_as_path}.json`;
    fs.mkdirSync(path.dirname(full_path), {recursive: true})
    fs.writeFileSync(
        full_path,
        JSON.stringify(schema, null, 4),
        "utf8",
    );
})
wrappedExamples.forEach((example) => {
    id_as_path = example["path"].replace("-", "_");
    full_path = `./docs/js/example/${id_as_path}.json`;
    fs.mkdirSync(path.dirname(full_path), {recursive: true})
    fs.writeFileSync(
        full_path,
        JSON.stringify(example["data"], null, 4),
        "utf8",
    );
})
