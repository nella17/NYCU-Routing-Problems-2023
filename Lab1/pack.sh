#!/bin/bash
rv () {
  echo -n "$1 : " >&2
  read x
  echo $x
}
id=`rv id`
rm -fr $id $id.zip
mkdir $id
cp ./*.cpp ./*.hpp Makefile $id
zip -r $id.zip $id
rm -fr $id
