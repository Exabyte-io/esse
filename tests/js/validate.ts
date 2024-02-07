/* eslint-disable no-unused-expressions */
// eslint-disable-next-line import/no-extraneous-dependencies
import { expect } from "chai";
import fs from "fs";
import groupBy from "lodash/groupBy";
import path from "path";

import { JSONSchemasInterface } from "../../src/js/esse/JSONSchemasInterface";
import { walkDirSync } from "../../src/js/scripts/utils";

const examplesPath = path.resolve("./lib/js/example");
const schemasPath = path.resolve("./lib/js/schema");

describe("validate all examples", () => {
    console.log("======");
    walkDirSync(examplesPath, (examplePath) => {
        if (examplePath.endsWith(".json")) {
            const schemaPath = examplePath.replace("/example/", "/schema/");
            const example = JSON.parse(fs.readFileSync(examplePath).toString());
            const schema = JSON.parse(fs.readFileSync(schemaPath).toString());

            console.log(schemaPath);

            const result = JSONSchemasInterface.validate(example, schema);

            if (!result.isValid) {
                console.log({
                    examplePath,
                    schemaPath,
                    errors: JSON.stringify(result.errors),
                });
            }

            expect(result.isValid).to.be.true;
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

    console.log(repeatedSchemaTitles);

    expect(repeatedSchemaTitles).to.be.an("array").that.is.empty;
});
