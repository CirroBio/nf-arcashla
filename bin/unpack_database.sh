#!/bin/bash

set -e

WORK_DIR=$(pwd)
REF_DIR=$(realpath /usr/local/share/arcas-hla-*/)
echo "Unpacking database to ${REF_DIR}/dat/IMGTHLA/"
mkdir -p ${REF_DIR}/dat/IMGTHLA/

tar -xzvf ${1}
mv IMGTHLA-*-alpha/* ${REF_DIR}/dat/IMGTHLA/
rm ${1}

cd ${REF_DIR}/dat/IMGTHLA/
for zip in *.zip; do
    unzip -o ${zip}
    rm ${zip}
done

cd ${WORK_DIR}
echo "\$(date) - Rebuilding database"
arcasHLA reference --rebuild --verbose 2>&1

echo "\$(date) - Database rebuilt"
