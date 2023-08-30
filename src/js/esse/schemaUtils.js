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

export function makeFlatSchemaKey(schemaId) {
    return schemaId.replace(/[^a-zA-Z-]/g, ":");
}

export function makeFlatSchemaRef(schemaId) {
    return { $ref: `#/definitions/${makeFlatSchemaKey(schemaId)}` };
}

export function buildSchemaDefinitions(originalSchemas) {
    const schemas = originalSchemas.map((schema) => {
        return mapObjectDeep(schema, (value) => {
            if (typeof value === "object" && value.$id) {
                return makeFlatSchemaRef(value.$id);
            }
            if (typeof value === "object" && value.$ref) {
                // assume value.$ref is a $id
                return makeFlatSchemaRef(value.$ref);
            }
        });
    });

    return keyBy(schemas, ({ $id }) => makeFlatSchemaKey($id));
}
