# exabyte-materials-json

Contains schemas and examples for materials and simulations related data in JSON representation. Also this repo represent npm 
package that can be used by clients in testing, on server and in some limited browser environments.

## Installation

```bash
cd <this repo dir>
pip install -r requirements.txt
```

To ensure latest functionality, run `pip install --upgrade https://github.com/timurbazhirov/json_include/archive/master.zip`.

## Usage

To produce json files with no inclusion statements:

```bash
python compile.py
```

`-m` flag will minify json files.

# Tests

Run from root directory of this repository:

```bash
npm install
npm test
```
