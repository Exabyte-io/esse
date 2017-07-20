import chai from "chai";
import Ajv from "ajv";
import path from "path";
import url from "url";
import _ from "lodash";
import prettyjson from "prettyjson";
import {NAMESPACE} from "../src/settings";
import {rawSchemas, examples, getSchemaByIdBasename, getExampleByIdBasename} from "../src/index";

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

describe('example', function () {
    examples.forEach(example => {
        it(`${example.id} should be valid`, function () {
            const validator = ajv.getSchema(example.id);
            const valid = validator(example);
            if (!valid) {
                console.log(prettyjson.render(validator.errors));
            }
            chai.expect(valid).to.be.ok;
        });
    });
});
//
//describe('unit', function () {
//    schemasUnit.forEach(schema => {
//        if (schema.example) {
//            const fullPath = path.join(schema.dir, schema.filename);
//            it(`${fullPath} should be valid`, function () {
//                const validator = ajvn.getSchema(schema.id);
//                const valid = validator(schema.example);
//                if (!valid) {
//                    console.log(prettyjson.render(validator.errors));
//                }
//                chai.expect(valid).to.be.ok;
//            });
//        }
//    });
//});
//
//describe('templates', function () {
//    templateList.forEach(tmpl => {
//        if (tmpl.template) {
//            it(`${tmpl.type} assignment should be valid`, function () {
//                const validator = ajv.getSchema("https://exabyte.io/schemas/job/model/method/workflow/unit.json");
//                const valid = validator(tmpl.template);
//                if (!valid) {
//                    console.log(prettyjson.render(validator.errors));
//                }
//                chai.expect(valid).to.be.ok;
//            });
//        }
//    });
//});
//
//describe('included unit schema', function () {
//    var ajvLocal = new Ajv({allErrors: true});
//    var validator = ajvLocal.compile(jsonInclude('job/model/method/workflow/unit.json'));
//    templateList.forEach(tmpl => {
//        if (tmpl.template) {
//            it(`${tmpl.type} assignment should be valid`, function () {
//                const valid = validator(tmpl.template);
//                if (!valid) {
//                    console.log(prettyjson.render(validator.errors));
//                }
//                chai.expect(valid).to.be.ok;
//            });
//        }
//    });
//});
