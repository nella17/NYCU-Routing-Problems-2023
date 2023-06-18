#!/bin/bash
if [ ! "$ID" ]; then
  echo "no id"
  exit
fi
rm -fr $ID.zip
zip -r $ID.zip Makefile src
