/* eslint-disable no-await-in-loop */
/* eslint-disable no-restricted-syntax */
import fsPromises from "fs/promises";
import fs from "fs";
import path from "path";

/**
 *
 * @param dir {String}
 * @param callback {() => void | Promise<void>}
 */
export async function walkDir(dir, callback) {
    const subDirs = await fsPromises.readdir(dir);

    for (const subDir of subDirs) {
        const itemPath = path.join(dir, subDir);
        const stat = await fsPromises.stat(itemPath);

        if (stat.isDirectory()) {
            await walkDir(itemPath, callback);
        } else {
            await callback(itemPath);
        }
    }
}

/**
 *
 * @param dir {String}
 * @param callback {() => void}
 */
export function walkDirSync(dir, callback) {
    const subDirs = fs.readdirSync(dir);

    for (const subDir of subDirs) {
        const itemPath = path.join(dir, subDir);
        const stat = fs.statSync(itemPath);

        if (stat.isDirectory()) {
            walkDirSync(itemPath, callback);
        } else {
            callback(itemPath);
        }
    }
}
