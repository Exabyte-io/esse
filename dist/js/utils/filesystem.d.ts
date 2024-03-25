export declare function walkDir(dir: string, callback: (itemPath: string) => void | Promise<void>): Promise<void>;
export declare function walkDirSync(dir: string, callback: (itemPath: string) => void): void;
