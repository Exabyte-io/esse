#!/usr/bin/env bash
set -e

PYTHON_BIN="python3"
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
VENV_NAME="venv"

usage() {
    echo "run-tests.sh [--disable-npm-install] [--venvdirv=venv] [--python-bin=/usr/bin/python3]"
    exit 1
}

check_args() {
    DISABLE_NPM_INSTALL="FALSE"
    for i in "$@"; do
        case $i in
            --disable-npm-install)
                DISABLE_NPM_INSTALL="TRUE"
                shift
                ;;
            -p=* | --python-bin=*)
                PYTHON_BIN="${i#*=}"
                ;;
            -v=* | --venvdir=*)
                VENV_NAME="${i#*=}"
                ;;
            *)
                usage
                ;;
        esac
    done
}

check_args $@

# Prepare the execution virtualenv
virtualenv --python ${PYTHON_BIN} ${THIS_DIR}/${VENV_NAME}
source ${THIS_DIR}/${VENV_NAME}/bin/activate
trap "deactivate" EXIT
if [ -f ${THIS_DIR}/requirements-dev.txt ]; then
    pip install -r ${THIS_DIR}/requirements-dev.txt --no-deps
fi

# python tests
coverage run -m esse.tests.validate

# Generate the code coverage reports
coverage report
coverage html --directory htmlcov_${TEST_TYPE}
coverage xml -o coverage_${TEST_TYPE}.xml

# Prepare the javascript environment
export NVM_DIR="/root/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && . "$NVM_DIR/nvm.sh"
if [ ${DISABLE_NPM_INSTALL} == "FALSE" ]; then
    npm install
fi

# Javascript tests
npm test
if [ $? -ne 0 ]; then
    echo "NPM TESTS FAILED!"
    exit 1
fi
