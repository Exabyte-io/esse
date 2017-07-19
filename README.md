# exabyte-materials-json

Contains schemas and examples for materials and simulations related data in JSON representation. Can be used as a node
or python module on server-side.

## Installation

Python:

```bash
cd <this repo dir>
virtualenv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

To ensure latest functionality, run `pip install --upgrade https://github.com/timurbazhirov/json_include/archive/master.zip`.

Node:
```bash
npm install
```

## Usage

To produce json files with no inclusion statements (python):

```bash
python compile.py
```

`-m` flag will minify json files.

To produce json array with all dereferenced schemas (javascript):

```bash
npm install
DEBUG=true node src/index.js > schemas.json
```

# Tests

Run from root directory of this repository:

```bash
npm test
```
