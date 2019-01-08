#!/usr/bin/env bash

SOURCE="${BASH_SOURCE[0]}"
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

if [ ! -d ${DIR}/.virtualenv ]; then
    virtualenv ${DIR}/.virtualenv
fi
source ${DIR}/.virtualenv/bin/activate
pip -q install -r ${DIR}/requirements.txt
export PYTHONPATH=${DIR}:${PYTHONPATH}

# javascript tests
npm install
npm test

# python tests
python tests/validate.py
