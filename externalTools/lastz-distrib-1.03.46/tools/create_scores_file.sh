#!/bin/bash
#
# create_scores_ path_to_encode_directories comparison_species

ENCODE="$1"
REGION="ENm010"
REFSPECIES="human"
SECSPECIES=$2

THISDIR=`dirname $0`

lastz_D --inferonly=${THISDIR}/create_scores_file.control \
    ${ENCODE}/${REGION}/${REFSPECIES}.${REGION}.fa \
    ${ENCODE}/${REGION}/${SECSPECIES}.${REGION}.fa \
  | ${THISDIR}/expand_scores_file.py --overridegaps

