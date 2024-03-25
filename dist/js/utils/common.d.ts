export declare function isInstanceOf(object: object, type: "Array"): object is object[];
export declare function isInstanceOf(object: object, type: "Object"): object is object;
/**
 * Makes sure that text representing arrays is parsed into Arrays.
 */
export declare function safeParseJSON(string: string): any;
