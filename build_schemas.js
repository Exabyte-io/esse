/**
 * build_schemas uses node API to read all schemas from the FS
 * at build time and appends them to the main index.js file for
 * downstream consumption to avoid FS calls in the browser.
 */
const fs = require("fs");
const { ESSE } = require("./lib/js/esse");

const { schemas } = new ESSE();

/**
 * Add additional export of all ESSE schemas at the end of the index.js
 * file to be able to import them in the same way as the ESSE class:
 * const { ESSE, schemas } from '@exabyte-io/esse.js';
 */
fs.appendFileSync(
    "./lib/js/esse/index.js",
    "\n\nexports.schemas = "+JSON.stringify(schemas)+";",
    "utf8",
);
