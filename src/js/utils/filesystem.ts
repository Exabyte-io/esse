/* eslint-disable no-await-in-loop */
/* eslint-disable no-restricted-syntax */
import fs from "fs";
import fsPromises from "fs/promises";
import path from "path";

export async function walkDir(dir: string, callback: (itemPath: string) => void | Promise<void>) {
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

export function walkDirSync(dir: string, callback: (itemPath: string) => void) {
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
