#!/bin/bash

source geometry.sh

echo "%%%%%%%%%%% STEP B04 - standard analysis of results, RDF %%%%%%%%%%%"

if [ ! -e solution.cfg ] ; then
  echo
  echo "ERROR! File solution.cfg is missing, repeat Step B01.sh"
  sleep 3
  exit
fi

res.sh solution rho
