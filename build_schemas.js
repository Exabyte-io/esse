/**
 * build_schemas uses node API to read all schemas from the FS
 * at build time and writes them out to a single schemas.js file for
 * downstream consumption to avoid FS calls in the browser.
 */
const fs = require("fs");
const { ESSE } = require("./lib/js/esse");

const { schemas } = new ESSE();
fs.appendFileSync(
    "./lib/js/esse/index.js",
    "\n\nvar _schemas = "+JSON.stringify(schemas)+"; \n\nexports.schemas = _schemas;",
    "utf8",
);
