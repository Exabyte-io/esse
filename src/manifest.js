import fs from "fs-extra";
import path from "path";
import {MANIFEST_DIR} from "./settings";
import YAML from "yamljs";
import _ from "lodash";

/**
 * This file contains functions for converting manifest data from yaml files to tree format:
 *       {
 *            "name" : "statistical",
 *            "slug" : "st",
 *            "children" : [
 *                {
 *                    "name" : "deterministic",
 *                    "slug" : "det",
 *                    "field" : "type",
 *                    "children" : [
 *                        {
 *                            "name" : "machine learning",
 *                            "slug" : "ml",
 *                            "field" : "subtype",
 *                            "children" : [
 *                                {
 *                                    "name" : "classification",
 *                                    "omit" : true,
 *                                    "slug" : "cl",
 *                                    "schema" : {
 *                                        "type" : "object"
 *                                    }
 *                                },
 *                                {
 *                                    "name" : "regression",
 *                                    "slug" : "re",
 *                                    "children" : null,
 *                                    "schema" : {
 *                                        "type" : "object"
 *                                    }
 *                                }
 *                            ]
 *                        }
 *                    ]
 *                }
 *            ]
 *        }
 */


/**
 * @summary Converts characteristic from schemas.yaml to manifest tree branch.
 * @param name {String}
 * @param ch {Object}
 * @return {Array}
 */
function characteristicToChildren(name, ch) {
    const children = [];
    if (ch.types) {
        ch.types.forEach(function (type) {
            children.push({
                slug: `${name}:${type}`,
                name: `${name}:${type}`,
                type: 'number'
            })
        });
    } else {
        children.push({
            slug: name,
            name: name,
            type: 'number'
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

/**
 * @summary Puts models array into tree structure.
 * @param props
 */
function models(props) {
    return {
        children: props
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
