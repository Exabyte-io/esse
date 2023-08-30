/* eslint-disable no-unused-expressions */
// eslint-disable-next-line import/no-extraneous-dependencies
import { expect } from "chai";
import file from "file";
import groupBy from "lodash/groupBy";
import path from "path";

import { ESSE } from "../index";
import { EXAMPLES_DIR, SCHEMAS_DIR } from "../settings";
import { parseIncludeReferenceStatements } from "../utils";

const esse = new ESSE();

describe("validate all examples", () => {
    file.walkSync(EXAMPLES_DIR, (dirPath_, dirs_, files_) => {
        files_.forEach((file_) => {
            const examplePath = path.join(dirPath_, file_);
            const schemaPath = examplePath.replace(EXAMPLES_DIR, SCHEMAS_DIR);
            if (!examplePath.endsWith(".json") || !schemaPath.endsWith(".json")) {
                // ignore files like .DS_Store
                return;
            }
            it(`${examplePath.replace(`${EXAMPLES_DIR}/`, "")}`, () => {
                const example = parseIncludeReferenceStatements(examplePath);
                const schema = parseIncludeReferenceStatements(schemaPath);
                const valid = esse.validate(example, schema);
                // eslint-disable-next-line no-unused-expressions
                expect(valid).to.be.ok;
            });
        });
    });
});

describe("schema titles must be unique or empty", () => {
    const repeatedSchemaTitles = Object.entries(groupBy(esse.schemas, "title"))
        .filter(([title, groupedValues]) => title !== "undefined" && groupedValues.length > 1)
        .map(([title, groupedValues]) => [title, groupedValues.map(({ $id }) => $id)]);

    console.log(repeatedSchemaTitles);

    expect(repeatedSchemaTitles).to.be.an("array").that.is.empty;
});
