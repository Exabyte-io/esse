import keyBy from "lodash/keyBy";

export function mapObjectDeep(object, mapValue) {
    if (typeof object !== "object") {
        return object;
    }

    if (Array.isArray(object)) {
        return object.map(
            (innerValue) => mapValue(innerValue) || mapObjectDeep(innerValue, mapValue),
        );
    }

    const entries = Object.entries(object).map(([key, value]) => {
        return [key, mapValue(value) || mapObjectDeep(value, mapValue)];
    });

    return Object.fromEntries(entries);
}

export function makeFlatSchemaId(schemaId) {
    return schemaId.replace(/[^a-zA-Z]/g, "-");
}

export function makeFlatSchemaRef(schemaId) {
    return { $ref: `#/definitions/${makeFlatSchemaId(schemaId)}` };
}

export function buildSchemaDefinitions(originalSchemas) {
    const schemas = originalSchemas.map((schema) => {
        return mapObjectDeep(schema, (value) => {
            if (typeof value === "object" && value.schemaId) {
                return makeFlatSchemaRef(value.schemaId);
            }
            if (typeof value === "object" && value.$ref) {
                // assume value.$ref is a schemaId
                return makeFlatSchemaRef(value.$ref);
            }
        });
    });

    return keyBy(schemas, ({ schemaId }) => makeFlatSchemaId(schemaId));
}
