"use strict";
var __importDefault =
    (this && this.__importDefault) ||
    function (mod) {
        return mod && mod.__esModule ? mod : { default: mod };
    };
Object.defineProperty(exports, "__esModule", { value: true });
exports.PROPERTIES_MANIFEST_PATH = exports.EXAMPLES_DIR = exports.SCHEMAS_DIR = void 0;
const path_1 = __importDefault(require("path"));
exports.SCHEMAS_DIR = path_1.default.resolve(__dirname, "../../../schema");
exports.EXAMPLES_DIR = path_1.default.resolve(__dirname, "../../../example");
exports.PROPERTIES_MANIFEST_PATH = path_1.default.resolve(
    __dirname,
    "../../../manifest/properties.yaml",
);
