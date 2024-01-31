/* eslint-disable no-await-in-loop */
/* eslint-disable no-restricted-syntax */
import fs from "fs/promises";
import path from "path";

/**
 *
 * @param dir {String}
 * @param callback {() => void | Promise<void>}
 */
export async function walkDir(dir, callback) {
    const subDirs = await fs.readdir(dir);

    for (const subDir of subDirs) {
        const itemPath = path.join(dir, subDir);
        const stat = await fs.stat(itemPath);

        if (stat.isDirectory()) {
            await walkDir(itemPath, callback);
        } else {
            await callback(itemPath);
        }
    }
}

export function patchSchema(object, patchObjectCallback) {
    if (Array.isArray(object)) {
        return object.map((item) => patchSchema(item, patchObjectCallback));
    }

    if (typeof object !== "object") {
        return object;
    }

    const cleanObject = patchObjectCallback(object);

    const entries = Object.entries(cleanObject).map(([key, value]) => {
        return [key, patchSchema(value, patchObjectCallback)];
    });

    return Object.fromEntries(entries);
}
