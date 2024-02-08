import { expect } from "chai";
import fs from "fs";
import groupBy from "lodash/groupBy";
import path from "path";

import * as ajv from "../../src/js/ajv/utils";
import { JSONSchemasInterface } from "../../src/js/esse/JSONSchemasInterface";
import { walkDirSync } from "../../src/js/scripts/utils";

const examplesPath = path.resolve("./lib/js/example");
const schemasPath = path.resolve("./lib/js/schema");

describe("validate all examples", () => {
    walkDirSync(examplesPath, (examplePath) => {
        if (examplePath.endsWith(".json")) {
            const schemaPath = examplePath.replace("/example/", "/schema/");
            const example = JSON.parse(fs.readFileSync(examplePath).toString());
            const schema = JSON.parse(fs.readFileSync(schemaPath).toString());

            console.log(`Validating example: ${examplePath}`);

            const result = ajv.validate(example, schema);

            if (!result.isValid) {
                console.log({
                    examplePath,
                    schemaPath,
                    errors: JSON.stringify(result.errors),
                    finalJSON: JSON.stringify(example),
                });
            }

            expect(result.isValid).to.be.equal(true);
        }
    });
});

describe("schema titles must be unique or empty", () => {
    JSONSchemasInterface.setSchemaFolder(schemasPath);
    const schemas = JSONSchemasInterface.schemasCache.values();
    const repeatedSchemaTitles = Object.entries(groupBy(schemas, "title"))
        .filter(([title, groupedValues]) => title !== "undefined" && groupedValues.length > 1)
        // @ts-ignore
        .map(([title, groupedValues]) => [title, groupedValues.map(({ $id }) => $id)]);

    if (repeatedSchemaTitles.length) {
        console.log(repeatedSchemaTitles);
    }

    // eslint-disable-next-line no-unused-expressions
    expect(repeatedSchemaTitles).to.be.an("array").that.is.empty;
});
