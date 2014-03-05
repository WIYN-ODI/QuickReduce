reset
echo
echo "Starting filterdef ..."
echo
../admin/filterdef.py odi_g.txt odi_r.txt odi_i.txt odi_z.txt quotaU.dat +wiynodi/lot6.txt +wiynodi/WIYN_primary.dat +wiynodi/WIYN_primary.dat +wiynodi/WIYN_primary.dat +wiynodi/ODI_coatings.txt +wiynodi/ODI_FusedSilica.txt +wiynodi/ODI_PBL6Y.txt +wiynodi/kpnotellurix.txt mosaic/*.txt +wiynodi/podifudge.txt +wiynodi/ADCSylgAndCoatingTrans.txt +wiynodi/ADCSylgAndCoatingTrans.txt | tee filterdefs

