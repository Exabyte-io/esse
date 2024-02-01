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
