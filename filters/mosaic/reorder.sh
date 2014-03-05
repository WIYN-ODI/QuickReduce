#!/bin/tcsh


sort -n $1 |sed '/^$/d' > temp
mv temp $1


