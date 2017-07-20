import Ajv from "ajv";
import chai from "chai";
import lodash from "lodash";
import prettyjson from "prettyjson";
import {rawSchemas, examples, getSchemaByIdBasename} from "../src/index";

const ajv = Ajv({
    allErrors: true,
    schemas: rawSchemas
});

function getSchemaAndExpectOK(basename) {
    const schema = getSchemaByIdBasename(basename);
    chai.expect(ajv.getSchema(schema.id)).to.be.ok;
}

describe('ajv', function () {
    it('works on primitive schema', function () {
        getSchemaAndExpectOK("1d_data_series.json")
    });
    it('works on schema with external ref', function () {
        getSchemaAndExpectOK("atomic_coordinates.json")
    });
});

describe('all examples are valid', function () {
    examples.forEach(example => {
        it(`${example.id} should be valid`, function () {
            const validator = ajv.getSchema(example.id);
            const valid = validator(example.content);
            if (!valid) {
                console.log(prettyjson.render(validator.errors));
            }
            chai.expect(valid).to.be.ok;
        });
    });
});
