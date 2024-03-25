"use strict";
var __importDefault =
    (this && this.__importDefault) ||
    function (mod) {
        return mod && mod.__esModule ? mod : { default: mod };
    };
Object.defineProperty(exports, "__esModule", { value: true });
exports.parseIncludeReferenceStatementsByDir = exports.parseIncludeReferenceStatements = void 0;
// @ts-ignore
const json_schema_deref_sync_1 = __importDefault(require("json-schema-deref-sync"));
const path_1 = __importDefault(require("path"));
const json_include_1 = require("../json_include");
const filesystem_1 = require("../utils/filesystem");
/**
 * Resolves `include` and `$ref` statements.
 * @param filePath {String} file to parse.
 */
function parseIncludeReferenceStatements(filePath) {
    const jsonResolver = new json_include_1.JSONInclude();
    const parsed = jsonResolver.parseIncludeStatements(filePath);
    const dirPath = path_1.default.dirname(filePath);
    let dereferenced = (0, json_schema_deref_sync_1.default)(parsed, {
        baseFolder: dirPath,
        removeIds: true,
    });
    // handle circular references and use non-dereferenced source
    if (dereferenced instanceof Error && dereferenced.message === "Circular self reference") {
        dereferenced = parsed;
    }
    return dereferenced;
}
exports.parseIncludeReferenceStatements = parseIncludeReferenceStatements;
/**
 * Resolves `include` and `$ref` statements for all the JSON files inside a given directory.
 * @param dirPath directory to parse.
 */
function parseIncludeReferenceStatementsByDir(dirPath, wrapInDataAndPath = false) {
    const schemas = [];
    const schemasWithPath = [];
    const topDir = path_1.default.resolve(__dirname, "../../../");
    (0, filesystem_1.walkDirSync)(dirPath, (filePath) => {
        if (filePath.endsWith(".json")) {
            const config = parseIncludeReferenceStatements(filePath);
            if (wrapInDataAndPath) {
                const _path = path_1.default.join(
                    // remove leading slashes and "example" from path
                    path_1.default
                        .dirname(filePath)
                        .replace(path_1.default.join(topDir, "example"), "")
                        .replace(/^\/+/, ""),
                    path_1.default.basename(filePath).replace(".json", ""),
                );
                schemasWithPath.push({ data: config, path: _path });
            } else {
                schemas.push(config);
            }
        }
    });
    return wrapInDataAndPath ? schemasWithPath : schemas;
}
exports.parseIncludeReferenceStatementsByDir = parseIncludeReferenceStatementsByDir;
