import { JSONSchema } from "./utils";

export type JSONSchemasInterfaceQuery = { [key in keyof JSONSchema]: { $regex: string } };

export default class JSONSchemasInterface {
    static schemasCache = new Map<string, JSONSchema>();

    static setSchemas(schema: JSONSchema[]) {
        schema.forEach((schema) => this.addSchema(schema));
    }

    static addSchema(schema: JSONSchema) {
        if (schema.$id) {
            this.schemasCache.set(schema.$id, schema);
        }
    }

    static getSchemaById(schemaId: string) {
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
    static matchSchema(query: JSONSchemasInterfaceQuery) {
        const searchFields = Object.keys(query) as Array<keyof typeof query>;

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
