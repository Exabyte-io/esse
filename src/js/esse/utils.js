import file from "file";
import deref from "json-schema-deref-sync";
import path from "path";

import { JSONInclude } from "../json_include";

/**
 * Resolves `include` and `$ref` statements.
 * @param filePath {String} file to parse.
 */
export function parseIncludeReferenceStatements(filePath) {
    const jsonResolver = new JSONInclude();
    const parsed = jsonResolver.parseIncludeStatements(filePath);
    const dirPath = path.dirname(filePath);
    let dereferenced = deref(parsed, { baseFolder: dirPath });
    // handle circular references and use non-dereferenced source
    if (dereferenced instanceof Error && dereferenced.message === "Circular self reference") {
        dereferenced = parsed;
    }
    return dereferenced;
}

/**
 * Resolves `include` and `$ref` statements for all the JSON files inside a given directory.
 * @param dirPath {String} directory to parse.
 */
export function parseIncludeReferenceStatementsByDir(dirPath, wrapInDataAndPath = false) {
    const data = [];
    const topDir = path.resolve(__dirname, "../../../");
    file.walkSync(dirPath, (dirPath_, dirs_, files_) => {
        files_.forEach((file_) => {
            const filePath = path.join(dirPath_, file_);
            if (filePath.endsWith(".json")) {
                const config = parseIncludeReferenceStatements(filePath);
                if (wrapInDataAndPath) {
                    const _path = path.join(
                        // remove leading slashes and "example" from path
                        path
                            .dirname(filePath)
                            .replace(path.join(topDir, "example"), "")
                            .replace(/^\/+/, ""),
                        path.basename(filePath).replace(".json", ""),
                    );
                    data.push({ data: config, path: _path });
                } else {
                    data.push(config);
                }
            }
        });
    });
    return data;
}
