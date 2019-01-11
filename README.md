# Exabyte Source of Schemas and Examples (ESSE)

ESSE contains schemas and examples for materials and simulations related data in JSON representation. 

## Installation

ESSE can be used as a Node or Python package on the server side.
Pleas note that schemas and examples are not available on the client.

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

ESSE can be install as a Node package either via NPM or the repository as below.

#### NPM

```bash
npm install exabyte-esse
```

#### Repository

Add `"exabyte-esse": "file:PATH_TO_ESSE_REPOSITIRY"` to `package.json`.

## Usage

ESSE contains separate but equivalent interfaces for Python and Node.
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

## Tests

Run the following command from the root directory of this repository to run the tests.
The script will run both Node and Python tests in which examples are validated against their schemas.

```bash
sh run-tests.sh
```

## Contribution

This is an open-source repository and we welcome contributions for other test cases.
We suggest forking this repository and introducing the adjustments there.
The changes in the fork can further be considered for merging into this repository as it is commonly used on Github.
This process is explained in more details [here](https://gist.github.com/Chaser324/ce0505fbed06b947d962).
