#!/bin/bash

# usage: ./makeCircosConf.sh <PAHT TO KARYOTYPE_FILE> <PATH TO INPUT PAIRS_DUMP> <PATH TO OUTPUT>
KARYOTYPE_FILE=$1  # First argument
LINK_FILE=$2       # Second argument
NEW_CONF=$3

# Use sed to replace placeholders
sed -i "s|__KARYOTYPE__|$KARYOTYPE_FILE|g" $NEW_CONF
sed -i "s|__LINK_FILE__|$LINK_FILE|g" $NEW_CONF