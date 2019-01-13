# ESSE

Exabyte Source of Schemas and Examples (ESSE) contains data formats and associated examples specifically designed for digital materials science [1](#links).

## Installation

ESSE can be used as a Node.js or Python package on the server side.
Please note that schemas and examples are not available on the client.

### Python

ESSE can be install as a Python package either via PyPi or the repository as below.

#### PyPi

```bash
pip install esse
```

#### Repository

```bash
virtualenv .venv
source .venv/bin/activate
pip install -e PATH_TO_ESSE_REPOSITIRY
```

### Node

ESSE can be install as a Node.js package either via NPM or the repository as below.

#### NPM

```bash
npm install exabyte-esse
```

#### Repository

Add `"exabyte-esse": "file:PATH_TO_ESSE_REPOSITIRY"` to `package.json`.

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
import {ESSE} from "exabyte-esse";

const es = new ESSE();
const schema = es.getSchemaById("material");
```

## Structure

ESSE contains 3 main directories, [schema](schema), [example](example) and [src](src) outlined below.

### Schema

The schema directory contains the schemas specifying the rules to structure materials-related data.
In order to apply object-oriented design principals, a set of core schemas, outlined below are defined to facilitate the schema modularity.

#### Primitive

[Primitive](schema/core/primitive) directory contains a set of custom primitives that extends default standard primitive types allowed by schema, such as String and Number.
Primitives are solely defined by the default primitives and can not be re-constructed from each other.

#### Abstract

[Abstract](schema/core/abstract) directory contains unit-less schemas that are constructed from default and custom primitives.

#### Reusable

[Reusable](schema/core/reusable) directory contains the schemas that are widely used in other schemas to avoid duplication, constructed from the abstract and primitive schemas.

#### Reference

[Reference](schema/core/reference) directory contains the schemas defining the rules to structure the references.

### Example

This directory contains the examples formed according to the schemas and implements the same directory structure as the schema directory.

### Src

This directory contains Python and Javascript interfaces implementing the functionality to access and validate schemas and examples.

## Tests

Execute the following command from the root directory of this repository to run the tests.
The script will run both Javascript and Python tests in which examples are validated against their schemas.

```bash
sh run-tests.sh
```

## Contribution

We welcome contributions for other test cases.
We suggest forking this repository and introducing the adjustments there.
The changes in the fork can further be considered for merging into this repository as it is commonly used on Github [#links](2).

## Best Practices

- Use unique IDs for schemas. One can run `sh refactor.sh` to automatically set the IDs and reformat examples.

- Do not use circular references in the schemas, instead leave the type as object and add explanation to description.

## Links

1: [Data Convention for Digital Materials Science](https://www.overleaf.com/project/5c240af344c4383e719ff286)

2: [GitHub Standard Fork & Pull Request Workflow](https://gist.github.com/Chaser324/ce0505fbed06b947d962)
