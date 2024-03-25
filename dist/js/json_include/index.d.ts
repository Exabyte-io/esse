export interface CacheObject {
    [key: string]: CacheObject | CacheObject[];
}
export declare class JSONInclude {
    JSON_INCLUDE_CACHE: CacheObject;
    /**
     * Extracts file name(path) from the include statement.
     * @param value string to extract the file name from.
     */
    _getIncludeFileName: (value: unknown) => string | null | undefined;
    /**
     * Walks a nested object to resolve include statements.
     * @param obj Object to traverse.
     * @param dirpath directory from which `obj` is obtained. Include statements are relative to `obj` path.
     */
    _walkObjectToInclude(obj: CacheObject | CacheObject[], dirpath: string): void;
    /**
     * Resolves `include` statements.
     * @param filePath file to parse.
     */
    parseIncludeStatements(filePath: string): CacheObject[];
}
