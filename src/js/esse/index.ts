/* eslint-disable class-methods-use-this */
import fs from "fs";
import yaml from "js-yaml";
import mergeAllOf from "json-schema-merge-allof";
import path from "path";

import { EXAMPLES_DIR, PROPERTIES_MANIFEST_PATH, SCHEMAS_DIR } from "./settings";
import { JSONSchema, JSONSchemaWithPath, parseIncludeReferenceStatementsByDir } from "./utils";

const SCHEMAS = parseIncludeReferenceStatementsByDir(SCHEMAS_DIR);
const EXAMPLES = parseIncludeReferenceStatementsByDir(EXAMPLES_DIR, true);
const PROPERTIES_MANIFEST = yaml.load(
    fs.readFileSync(PROPERTIES_MANIFEST_PATH, { encoding: "utf-8" }),
) as object;
const RESULTS = Object.entries(PROPERTIES_MANIFEST)
    .map((k) => (k[1].isResult ? k[0] : null))
    .filter((x) => x) as object;

interface EsseConfig {
    schemas: JSONSchema[];
    examples: JSONSchema[];
    wrappedExamples: JSONSchemaWithPath[];
    propertiesManifest: object;
    results: object;
}

export class ESSE implements EsseConfig {
    readonly schemas: JSONSchema[];

    readonly examples: JSONSchema[];

    readonly wrappedExamples: JSONSchemaWithPath[];

    readonly propertiesManifest: object;

    readonly results: object;

    constructor(config?: EsseConfig) {
        this.schemas = config?.schemas || SCHEMAS;
        this.wrappedExamples = config?.wrappedExamples || EXAMPLES;
        this.examples = this.wrappedExamples.map((example) => example.data);
        this.propertiesManifest = config?.propertiesManifest || PROPERTIES_MANIFEST;
        this.results = config?.results || RESULTS;
    }

    writeResolvedSchemas(subfolder: string, skipMergeAllOff = false) {
        const schemasFolder = `${subfolder}/schema`;
        fs.rmSync(schemasFolder, { recursive: true, force: true });
        this.schemas.forEach((s) => {
            let mergedSchema = s;
            if (!skipMergeAllOff) {
                mergedSchema = mergeAllOf(s, {
                    resolvers: { defaultResolver: mergeAllOf.options.resolvers.title },
                });
            }
            const id_as_path = mergedSchema.$id?.replace(/-/g, "_");
            const full_path = `${schemasFolder}/${id_as_path}.json`;
            fs.mkdirSync(path.dirname(full_path), { recursive: true });
            fs.writeFileSync(full_path, JSON.stringify(mergedSchema, null, 4), "utf8");
        });
    }

    writeResolvedExamples(subfolder: string) {
        const examplesFolder = `${subfolder}/example`;
        fs.rmSync(`${examplesFolder}`, { recursive: true, force: true });
        this.wrappedExamples.forEach((e) => {
            const id_as_path = e.path.replace(/-/g, "_");
            const full_path = `${examplesFolder}/${id_as_path}.json`;
            fs.mkdirSync(path.dirname(full_path), { recursive: true });
            fs.writeFileSync(full_path, JSON.stringify(e.data, null, 4), "utf8");
        });
    }
}
