import fs from "fs-extra";
import path from "path";
import _ from "lodash";
import YAML from "yamljs";
import {MANIFEST_DIR} from "./settings";

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
 * @summary Converts characteristic from manifest representation to tree branch.
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
 * @summary Recursively scans nested object and omits all childred array elements with `omit` flag set to true
 * @param {Object} obj
 */
function _falseIfOmit(obj) {
    if (obj.omit) {
        return false;
    } else {
        if (obj.children) {
            obj.children = _.filter(obj.children, el => _falseIfOmit(el))
        }
    }
    return true;
}

/**
 * @summary Puts models array into tree structure.
 * @param props
 */
function models(props) {
    const modelsWithoutOmitted = _.filter(props, el => _falseIfOmit(el));
    return {
        children: modelsWithoutOmitted
    };
}

const MANIFESTS = [
    // `id` key is used further to retrieve the manifest data
    {
        id: 'properties',
        basename: 'properties.yaml',
        mapper: characteristics
    },
    {
        id: 'models',
        basename: 'models.yaml',
        mapper: models
    },
    {
        id: 'apps',
        basename: 'apps.yaml',
        mapper: x => x
    }
].map(manifest => {
    const filePath = path.join(MANIFEST_DIR, manifest.basename);
    const content = YAML.parse(fs.readFileSync(filePath).toString('utf8'));
    return Object.assign(
        manifest.mapper(content),
        {
            id: manifest.id,
            // relative path
            path: filePath.replace(`${MANIFEST_DIR}`, '').replace(/^\//, '')
        }
    );
});

export default MANIFESTS;
