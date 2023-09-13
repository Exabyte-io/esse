[![PyPI version](https://badge.fury.io/py/esse.svg)](https://badge.fury.io/py/esse)
[![npm version](https://badge.fury.io/js/%40exabyte-io%2Fesse.js.svg)](https://badge.fury.io/js/%40exabyte-io%2Fesse.js)
[![License: Apache](https://img.shields.io/badge/License-Apache-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

# ESSE

Essential Source of Schemas and Examples (ESSE) contains data formats and associated examples specifically designed for digital materials science (see refs. [1, 2](#links) below).

## Installation

ESSE can be used as a Node.js or Python package on the server side.

### Python

ESSE is compatible with Python 3.6+.  It can be installed as a Python package either via PyPI or the repository as below.

#### PyPI

```bash
pip install esse
```

#### Repository

```bash
virtualenv .venv
source .venv/bin/activate
pip install -e PATH_TO_ESSE_REPOSITORY
```

### Node

ESSE can be installed as a Node.js package either via NPM or the repository as below.

#### NPM

```bash
npm install @exabyte-io/esse.js
```

#### Repository

Add `"esse-js": "file:PATH_TO_ESSE_REPOSITORY"` to `package.json`.

## Usage

ESSE contains separate but equivalent interfaces for Python and Javascript.
The package provides `ESSE` class that can be initialized and used as below.

### Python

```python
from esse import ESSE

es = ESSE()
schema = es.get_schema_by_id("material")
```

### Node

```javascript
import {ESSE} from "esse-js";

const es = new ESSE();
const schema = es.getSchemaById("material");
```

## Structure

ESSE contains 3 main directories, [schema](schema), [example](example) and [src](src) outlined below.

### Schema

The schema directory contains the schemas specifying the rules to structure data. A set of core schemas, outlined below, are defined to facilitate the schema modularity.

#### Primitive

[Primitive](schema/core/primitive) directory contains a set of custom primitives that extends default standard primitive types allowed by schema, such as String and Number.
Primitives are solely defined by the default primitives and can not be re-constructed from each other.

#### Abstract

[Abstract](schema/core/abstract) directory contains unit-less schemas that are constructed from default and custom primitives.

#### Reusable

[Reusable](schema/core/reusable) directory contains the schemas that are widely used in other schemas to avoid duplication, constructed from the abstract and primitive schemas.

#### Reference

[Reference](schema/core/reference) directory contains the schemas defining the rules to structure the references to data sources.

### Example

This directory contains the examples formed according to the schemas and implements the same directory structure as the schema directory.

### src

This directory contains Python and Javascript interfaces implementing the functionality to access and validate schemas and examples.

### Generative vs Non-generative keys
Generative keys are the fields which allow for user input prior to calculation of the final property values. A flag is included in the schema comments on the fields in [property schemas](schema/properties_directory): `isGenerative:true` marks which fields to use as subschemas in the generation of a user input schema.
- On properties allowing user inputs, additional fields may be tagged, as in [the `file_content` property](schema/properties_directory/non-scalar/file_content.json)

## Tests

Execute the following command from the root directory of this repository to run the tests. The script will run both Javascript and Python tests in which examples are validated against the corresponding schemas.

```bash
bash run-tests.sh
```
The script has been tested with node.js v12.16.3 and v8.17.0 as well as Python version 2.7 (up to version 2.3.0) and 3.6+ (for version 2020.10.19 and later).

## Contribution

This repository is an [open-source](LICENSE.md) work-in-progress and we welcome contributions. We suggest forking this repository and introducing the adjustments there, the changes in the fork can further be considered for merging into this repository as it is commonly done on Github (see [3](#links) below).

## Best Practices

- Use unique IDs for schemas. One can run `sh refactor.sh` to automatically set the IDs and reformat examples.

- Do not use circular references in the schemas, instead leave the type as object and add explanation to description.

## Links

1: [Data-centric online ecosystem for digital materials science](https://arxiv.org/pdf/1902.10838.pdf)

2: [CateCom: A Practical Data-Centric Approach to Categorization of Computational Models](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.2c00112)

3: [GitHub Standard Fork & Pull Request Workflow](https://gist.github.com/Chaser324/ce0505fbed06b947d962)
