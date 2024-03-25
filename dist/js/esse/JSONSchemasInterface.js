"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
class JSONSchemasInterface {
    static setSchemas(schema) {
        schema.forEach((schema) => this.addSchema(schema));
    }
    static addSchema(schema) {
        if (schema.$id) {
            this.schemasCache.set(schema.$id, schema);
        }
    }
    static getSchemaById(schemaId) {
        return this.schemasCache.get(schemaId);
    }
    /**
     * @example <caption>Search by $id regex</caption>
     * JSONSchemasInterface.matchSchema({
     *   $id: {
     *     $regex: 'software-application'
     *   }
     * })
     *
     * @example <caption>Search by $id and title regex</caption>
     * JSONSchemasInterface.matchSchema({
     *   $id: {
     *     $regex: 'software-application'
     *   },
     *   title: {
     *     $regex: 'application'
     *   }
     * })
     */
    static matchSchema(query) {
        const searchFields = Object.keys(query);
        return Array.from(this.schemasCache.values()).find((schema) => {
            return searchFields.every((field) => {
                const queryField = query[field];
                const schemaField = schema[field];
                if (!queryField || typeof schemaField !== "string") {
                    return;
                }
                return new RegExp(queryField.$regex).test(schemaField);
            });
        });
    }
}
exports.default = JSONSchemasInterface;
JSONSchemasInterface.schemasCache = new Map();
