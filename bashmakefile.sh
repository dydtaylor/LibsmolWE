#!/bin/bash
mkdir -p ../wesmoldynresults
DATE=$(date +"%Y%m%d%H%M%S")
HASH=$(git rev-parse HEAD)
PREFIX="ReretestN5BrianUnitless"
cp -r mfiles "../wesmoldynresults/${PREFIX}_date_${DATE}hash_${HASH}"
