
dir1=/home/lyl/WORK3/qinyi/double-ITCZ/data/obs/
dir2=/home/lyl/WORK3/qinyi/double-ITCZ/data/atm/

case=(ISCCP BC5_f19g16 BC5_f19g16_mac2)
append=(_MONTHS .cam.h0.MONTHS .cam.h0.MONTHS)
dir=($dir1 $dir2 $dir2)

ncase=${#case[@]}

export lonS=230. #270.
export lonE=255. #290.
export latS=-20. #-30.
export latE=-5.  #-10.

#BC5_f19g16_mac2.cam.h0.MONTHS.nc
#BC5_f19g16.cam.h0.MONTHS.nc
#ISCCP_MONTHS.nc

workdir=/home/lyl/WORK3/qinyi/double-ITCZ/
cd $workdir


# form the regrid information
#tmp=CERES-EBAF_ANN_climo.nc
#cdo griddes $dir1/$tmp > regrid.txt
# note: here need to delete information of grid2 in regrid.txt

for icase in `seq 0 $[$ncase-1]`
do
	# extract cloud variables from dataset
	ncks -O -v CLDLOW ${dir[icase]}/${case[icase]}${append[icase]}.nc ${dir[icase]}/${case[icase]}${append[icase]}_CLDLOW.nc
	# rename to make them put into one file and not conflict with each other
	ncrename -v CLDLOW,CLDLOW$icase ${dir[icase]}/${case[icase]}${append[icase]}_CLDLOW.nc
	# regrid data
	cdo remapbil,regrid.txt ${dir[icase]}/${case[icase]}${append[icase]}_CLDLOW.nc ${dir[icase]}/${case[icase]}${append[icase]}_CLDLOW_regrid.nc
	# subsetting SEP region: (270-290,-30,-10)
	# notion: if you want to use real longitude value, please use "230." rather than "230".
	# if you write "230", it thinks it is the index rather than the real value.
	ncks -O -d lon,$lonS,$lonE -d lat,$latS,$latE ${dir[icase]}/${case[icase]}${append[icase]}_CLDLOW_regrid.nc ${dir[icase]}/${case[icase]}${append[icase]}_CLDLOW_regrid_SEP.nc
	# area average
	ncwa -O -a lat,lon ${dir[icase]}/${case[icase]}${append[icase]}_CLDLOW_regrid_SEP.nc ${dir[icase]}/${case[icase]}${append[icase]}_CLDLOW_regrid_SEP_avg.nc
	# stitch low cloud fraction from three cases into one file
	ncks -A ${dir[icase]}/${case[icase]}${append[icase]}_CLDLOW_regrid_SEP_avg.nc -o $workdir/CLDLOW_3cases_SEP_avg.nc
done

ncl season-value-lowcloud-SEP.ncl
