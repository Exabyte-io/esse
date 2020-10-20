#!/usr/bin/env bash

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"

if [ ! -d ${THIS_DIR}/venv ]; then
    virtualenv ${THIS_DIR}/venv
fi
source ${THIS_DIR}/venv/bin/activate
pip install -r ${THIS_DIR}/requirements-dev.txt --no-deps

python -c 'from esse.utils import *; refactor_examples(); refactor_schemas()'
