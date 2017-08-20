import fs from "fs-extra";
import path from "path";
import {MANIFEST_DIR} from "./settings";
import YAML from "yamljs";
import _ from "lodash";

function characteristicToChildren(name, ch) {
    const children = [];
    if (ch.types) {
        ch.types.forEach(function (type) {
            children.push({
                slug: `${name}:${type}`,
                name: `${name}:${type}`,
                type: 'number' // TODO: get from schema
            })
        });
    } else {
        children.push({
            slug: name,
            name: name,
            type: 'number' // TODO: get from schema
        });
    }

    return children;
}

function characteristics(props) {
    const result = [];

    _.each(props, function (value, name) {
        // return only characteristic properties
        if (value.characteristic) {
            result.push(...characteristicToChildren(name, value))
        }
    });

    return {
        children: result
    };
}

function models(props) {
    const result = props;
    return {
        children: result
    };
}

const MANIFEST_LIST = [
    {file: 'schemas.yaml', mapper: characteristics},
    {file: 'models.yaml', mapper: models}
];

export default MANIFEST_LIST.map(manifest=> {
    const filePath = path.join(MANIFEST_DIR, manifest.file);
    const content = YAML.parse(fs.readFileSync(filePath).toString('utf8'));
    return Object.assign(
        manifest.mapper(content),
        {id: manifest.file, path: filePath}
    );
});
