// eslint-disable-next-line import/no-extraneous-dependencies
import { expect } from "chai";
import file from "file";
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
