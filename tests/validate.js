import file from "file";
import path from "path";
import chai from "chai";

import {ESSE} from "../src/index";
import {EXAMPLES_DIR, SCHEMAS_DIR} from "../src/settings";
import {JSONSchemaResolver} from "../src/resolver/resolver";

const esse = new ESSE();
const jsonResolver = new JSONSchemaResolver();


describe('validate all examples', function () {
    file.walkSync(EXAMPLES_DIR, (dirPath_, dirs_, files_) => {
        files_.forEach(file_ => {
            const examplePath = path.join(dirPath_, file_);
            const schemaPath = examplePath.replace(EXAMPLES_DIR, SCHEMAS_DIR);
            it(`${examplePath.replace(`${EXAMPLES_DIR}/`, "")}`, function () {
                const example = jsonResolver.resolve(examplePath);
                const schema = jsonResolver.resolve(schemaPath);
                const valid = esse.validate(example, schema);
                chai.expect(valid).to.be.ok;
            })
        });
    });
});
