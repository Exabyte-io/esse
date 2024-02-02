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
        const res = mapValue(value);

        return [key, res === undefined ? mapObjectDeep(value, mapValue) : res];
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

export function walkSchema(object, callback) {
    if (Array.isArray(object)) {
        return object.map((item) => walkSchema(item, callback));
    }

    const cleanObject = callback(object);

    if (typeof cleanObject !== "object") {
        return cleanObject;
    }

    const entries = Object.entries(cleanObject).map(([key, value]) => {
        return [key, walkSchema(value, callback)];
    });

    return Object.fromEntries(entries);
}

export function addAdditionalPropertiesToSchema(schema, additionalProperties = false) {
    return walkSchema(schema, (object) => {
        if (typeof object !== "object") {
            return object;
        }

        if (object.type === "object" && object.properties && !("additionalProperties" in object)) {
            return {
                ...object,
                additionalProperties,
            };
        }
        return object;
    });
}
