#!/bin/sh
grep -v '#' odi_NB422.raw | tac | awk '{f=$2<=0?1e-5:$2; print $1*10, f; }' > odi_NB422b.txt
grep -v '#' odi_NB659.raw | tac | awk '{f=$2<=0?1e-5:$2; print $1*10, f; }' > odi_NB659.txt
grep -v '#' odi_NB746.raw | tac | awk '{f=$2<=0?1e-5:$2; print $1*10, f; }' > odi_NB746.txt

../admin/filterdef.py odi_NB746.txt odi_NB659.txt odi_NB422.txt odi_NB422b.txt  +wiynodi/lot6.txt +wiynodi/WIYN_primary.dat +wiynodi/WIYN_primary.dat +wiynodi/WIYN_primary.dat +wiynodi/ODI_coatings.txt +wiynodi/ODI_FusedSilica.txt +wiynodi/ODI_PBL6Y.txt +wiynodi/kpnotellurix.txt +wiynodi/podifudge.txt +wiynodi/ADCSylgAndCoatingTrans.txt +wiynodi/ADCSylgAndCoatingTrans.txt


