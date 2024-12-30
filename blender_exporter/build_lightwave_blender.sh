#!/bin/bash
set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
SRC_DIR=$(dirname ${SCRIPT_DIR})
CUR_DIR=$(pwd)

NAME=lightwave_blender

BOB_DIR=${SRC_DIR}/deps/bob
LIB_DIR=${SCRIPT_DIR}/${NAME}

cd ${BOB_DIR}
python3.11 setup.py build_ext --build-lib=${LIB_DIR}

cd ${SCRIPT_DIR}

zip -qr ${NAME} ${NAME}

echo "Blender plugin built. Open Blender and go to 'Edit - Preferences - Addons - Install...'"
echo "Install '${SCRIPT_DIR}/${NAME}.zip'."

cd ${CUR_DIR}