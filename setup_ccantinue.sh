#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export CCNUEROOT=$DIR
export PYTHONPATH=$CCNUEROOT:$PYTHONPATH
# Alternative configuration folder
if [ -z "$1" ]; then
config="CCAntiNuE"
else 
config=$1
fi
echo $config
export CCNuEConfig=$config
export PYTHONPATH=${CCNUEROOT}/configs/${config}:$PYTHONPATH
