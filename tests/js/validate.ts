import { expect } from "chai";
import fs from "fs";
import groupBy from "lodash/groupBy";
import path from "path";

import JSONSchemasInterface from "../../src/js/esse/JSONSchemasInterfaceServer";
import { AnyObject } from "../../src/js/esse/types";
import * as ajv from "../../src/js/utils/ajv";
import { walkDirSync } from "../../src/js/utils/filesystem";

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

interface Example extends AnyObject {
    property: string | Date;
}

describe("validate Date object", () => {
    const schema = {
        $id: "validate-date-object",
        type: "object",
        properties: {
            property: {
                type: "string",
                format: "date-time",
            },
        },
    };

    const example1: Example = { property: new Date() };
    const example2: Example = { property: "December 17, 1995 03:24:00" };
    const example3: Example = { property: "Invalid Date" };

    const result1 = ajv.validate(example1, schema);
    const result2 = ajv.validate(example2, schema);
    const result3 = ajv.validate(example3, schema);

    expect(result1.isValid).to.be.equal(true);
    expect(result2.isValid).to.be.equal(true);
    expect(result3.isValid).to.be.equal(false);

    expect(example1.property instanceof Date).to.be.equal(true);
    expect(example2.property instanceof Date).to.be.equal(true);
    expect(typeof example3.property === "string").to.be.equal(true);
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
