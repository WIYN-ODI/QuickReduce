#!/bin/bash

#
# This script tests the installed WCS library for a particular bug 
# that will cause the best-fit WCS solution to randomly be off by a 
# small amount (few arcseconds, mostly in declination)
#
# (c) 2014 Ralf Kotulla, All rights reserved
#

tries=25
echo "Running $tries WCS tests"
echo
echo "Correct solution is:"
echo "[ 202.46897694   47.13329625]"
echo "------------------------------"
echo
echo "Test results: ________"

if [ -e test.results ]
then
    # Delete the file test.results
    rm test.results
fi

#!/bin/bash
for (( c=1; c<=$tries; c++ ))
do
    # echo "Welcome $c times"
    printf "Try %02d: " $c
    ./test.py ota33.cat ota33.head | tee -a test.results
done

# for i in {1..$tries}
# do
#    echo "Welcome $i times"
# done

nresults=`cat test.results | uniq | wc -l`

echo
echo "Summary:"
echo "----------------------------------"
echo "Number of different results: $nresults"

if [ $nresults == 1 ]
then
   echo "==> Test successful, WCS works the way it should";
else
   echo "==> WCS has problems, use at your own risk!";
fi

if [ -e test.results ]
then
    # Delete the file test.results
    rm test.results
fi


echo
