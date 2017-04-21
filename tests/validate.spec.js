import chai from "chai";
import Ajv from "ajv";
import path from "path";
import url from "url";
import _ from "lodash";
import {getNormalizedSchemas} from "../src/json-schema-normalize";
import {templateList, schemasList, jsonInclude} from "../src/unit-schema-include";
import prettyjson from "prettyjson";

var NAMESPACE = "https://exabyte.io/schemas/";
var schemasUnit = schemasList;

// create the array of normalized schema descriptions
var schemas = getNormalizedSchemas(true);
var ajv = Ajv({
    allErrors: true,
    schemas: schemas.map(schema => schema.content)
});
var ajvn = Ajv({
    allErrors: true,
    schemas: schemasUnit.map(schema => schema.content)
});

schemasUnit.forEach((rawSchema) => {
    const sameUrlSchema = _.find(schemas, (normalizedSchema) => {
        var id = url.resolve(NAMESPACE, path.join(rawSchema.dir, rawSchema.filename));
        return normalizedSchema.id === id;
    });
    if (sameUrlSchema && sameUrlSchema.example) {
        rawSchema.example = sameUrlSchema.example;
    }
});

describe('ajv', function () {
    it('works on primitive schema', function () {
        chai.expect(ajv.getSchema('https://exabyte.io/schemas/primitive/1d_data_series.json')).to.be.ok;
    });

    it('works on schema with external ref', function () {
        chai.expect(ajv.getSchema('https://exabyte.io/schemas/material/structure/atomic_coordinates.json')).to.be.ok;
    });
});

describe('example', function () {
    schemas.forEach(schema => {
        if (schema.example) {
            const fullPath = path.join(schema.dir, schema.filename);
            it(`${fullPath} should be valid`, function () {
                const validator = ajv.getSchema(schema.id);
                const valid = validator(schema.example);
                if (!valid) {
                    console.log(JSON.stringify(schema.example))
                    console.log(prettyjson.render(validator.errors));
                }
                chai.expect(valid).to.be.ok;
            });
        }
    });
});

describe('unit', function () {
    schemasUnit.forEach(schema => {
        if (schema.example) {
            const fullPath = path.join(schema.dir, schema.filename);
            it(`${fullPath} should be valid`, function () {
                const validator = ajvn.getSchema(schema.id);
                const valid = validator(schema.example);
                if (!valid) {
                    console.log(prettyjson.render(validator.errors));
                }
                chai.expect(valid).to.be.ok;
            });
        }
    });
});

describe('templates', function () {
    templateList.forEach(tmpl => {
        if (tmpl.template) {
            it(`${tmpl.type} assignment should be valid`, function () {
                const validator = ajv.getSchema("https://exabyte.io/schemas/job/model/method/workflow/unit.json");
                const valid = validator(tmpl.template);
                if (!valid) {
                    console.log(prettyjson.render(validator.errors));
                }
                chai.expect(valid).to.be.ok;
            });
        }
    });
});

describe('included unit schema', function () {
    var ajvLocal = new Ajv({allErrors: true});
    var validator = ajvLocal.compile(jsonInclude('job/model/method/workflow/unit.json'));
    templateList.forEach(tmpl => {
        if (tmpl.template) {
            it(`${tmpl.type} assignment should be valid`, function () {
                const valid = validator(tmpl.template);
                if (!valid) {
                    console.log(prettyjson.render(validator.errors));
                }
                chai.expect(valid).to.be.ok;
            });
        }
    });
});
