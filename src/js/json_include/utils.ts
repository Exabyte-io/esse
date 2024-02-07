export function isInstanceOf(object: object, type: "Array"): object is object[];

export function isInstanceOf(object: object, type: "Object"): object is object;

export function isInstanceOf(object: object, type: string): boolean {
    return Object.prototype.toString.call(object).slice(8, -1) === type;
}

/**
 * Makes sure that text representing arrays is parsed into Arrays.
 */
export function safeParseJSON(string: string) {
    const obj = JSON.parse(string);
    return string[0] === "[" ? Object.keys(obj).map((key) => obj[key]) : obj;
}
