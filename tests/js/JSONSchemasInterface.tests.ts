import { assert, expect } from "chai";
import * as path from "path";

import allSchemas from "../../dist/js/schemas.json";
import JSONSchemasInterface from "../../src/js/esse/JSONSchemasInterfaceServer";
import { JSONSchema } from "../../src/js/esse/utils";

function assertSystemInSetSchema(schema?: JSONSchema) {
    const inSet = schema?.properties?.inSet as JSONSchema | undefined;
    const inSetItems = inSet?.items as JSONSchema | undefined;

    expect(schema).to.be.an("object");
    assert(schema?.$id, "system/in-set");
    expect(inSetItems?.properties?._id).to.be.an("object");
    expect(inSetItems?.properties?.cls).to.be.an("object");
    expect(inSetItems?.properties?.slug).to.be.an("object");
    expect(inSetItems?.properties?.type).to.be.an("object");
    expect(inSetItems?.properties?.index).to.be.an("object");
}

describe("JSONSchemasInterfaceServer", () => {
    it("can find schemas from esse dist folder; the schema is merged and clean", async () => {
        const schema = JSONSchemasInterface.getSchemaById("system/in-set");
        assertSystemInSetSchema(schema);
    });

    it("can find registered schemas; the schema is merged and clean", async () => {
        JSONSchemasInterface.setSchemaFolder(path.join(__dirname, "./fixtures/json"));

        const schema = JSONSchemasInterface.getSchemaById("system/in-set");
        assertSystemInSetSchema(schema);
    });
});

describe("JSONSchemasInterface", () => {
    it("can find registered schemas; the schema is merged and clean", async () => {
        JSONSchemasInterface.setSchemas(allSchemas as JSONSchema[]);

        const schema = JSONSchemasInterface.getSchemaById("system/in-set");
        assertSystemInSetSchema(schema);
    });
});
