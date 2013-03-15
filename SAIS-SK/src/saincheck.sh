#!/bin/sh

set -e -x

for filename in `ls *`
do
  sk-sain.x $filename
done
