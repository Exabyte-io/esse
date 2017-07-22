import chai from "chai";
import prettyjson from "prettyjson";
import {ajvHandler as ajv, examples, getSchemaByIdBasename} from "../src/index";

function getSchemaAndExpectOK(basename) {
    const schema = getSchemaByIdBasename(basename);
    chai.expect(ajv.getSchema(schema.id)).to.be.ok;
}

describe('test handpicked schemas', function () {
    it('primitive schema', function () {
        getSchemaAndExpectOK("1d_data_series.json")
    });
    it('schema with external ref', function () {
        getSchemaAndExpectOK("atomic_coordinates.json")
    });
});

describe('validate all examples', function () {
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
