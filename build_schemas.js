/**
 * build_schemas uses node API to read all schemas from the FS
 * at build time and writes them out to a single schemas.js file for
 * downstream consumption to avoid FS calls in the browser.
 */
 const fs = require("fs");
 const {ESSE} = require("./lib/js/esse");
 const {schemas} = new ESSE();
 fs.writeFileSync("./lib/js/esse/schemas.js", "module.exports = {schemas: " + JSON.stringify(schemas) + "}", "utf8");
