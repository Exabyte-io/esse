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
export NVM_DIR="/root/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && . "$NVM_DIR/nvm.sh"

npm install
npm test

if [ $? -ne 0 ]; then
    echo "NPM TESTS FAILED!"
    exit 1
fi

# python tests
python tests/validate.py

if [ $? -ne 0 ]; then
    echo "PYTHON TESTS FAILED!"
    exit 1
fi
