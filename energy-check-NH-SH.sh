
dir1=/home/lyl/WORK3/qinyi/double-ITCZ/data/obs/
dir2=/home/lyl/WORK3/qinyi/double-ITCZ/data/atm/

var1=(FSNS,FLNS,LHFLX,SHFLX,FSNT,FLUT)
var2=(FSNS,FLNS,LHFLX,SHFLX,FSNT,FLNT,PRECC,PRECL,PRECSC,PRECSL,QFLX)


case=(JRA25 BC5_f19g16 BC5_f19g16_mac2)
append=(_ANN_climo _ANN_climo _ANN_climo)
dir=($dir1 $dir2 $dir2)
var=($var1 $var2 $var2)

ncase=${#case[@]}

# SH
export lonS1=0.
export lonE1=360.
export latS1=-90.
export latE1=0.  

# NH
export lonS2=0.
export lonE2=360.
export latS2=0.
export latE2=90. 


workdir=/home/lyl/WORK3/qinyi/double-ITCZ/
cd $workdir


# form the regrid information
#tmp=CERES-EBAF_ANN_climo.nc
#cdo griddes $dir1/$tmp > regrid.txt
# note: here need to delete information of grid2 in regrid.txt

for icase in `seq 0 $[$ncase-1]`
do
	# extract cloud variables from dataset
	ncks -O -v ${var[icase]} ${dir[icase]}/${case[icase]}${append[icase]}.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeatvars.nc
	# calculate net heat flux over NH and SH
	if [ $icase == 0 ]; then
		ncap2 -O -s 'NetAtmHeat=((-1.0)*FSNT-FLUT)-((-1.)*FSNS-FLNS-LHFLX-SHFLX)' ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeatvars.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc
		ncap2 -O -s 'NetTOAHeat=((-1.0)*FSNT-FLUT)' ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc
		ncap2 -O -s 'NetSfcHeat=((-1.0)*FSNS-FLNS-LHFLX-SHFLX)' ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc
	else
		ncap2 -O -s 'LHFLX_new=(2501000.+333700)*QFLX-333700*1.e3*(PRECC+PRECL-PRECSC-PRECSL)' ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeatvars.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeatvars.nc
		ncap2 -O -s 'NetAtmHeat=(FSNT-FLNT)-(FSNS-FLNS-LHFLX_new-SHFLX)' ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeatvars.nc  ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc
		ncap2 -O -s 'NetTOAHeat=(FSNT-FLNT)' ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc  ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc
		ncap2 -O -s 'NetSfcHeat=(FSNS-FLNS-LHFLX_new-SHFLX)' ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc  ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc
	fi

	# regrid data
	cdo remapbil,regrid.txt ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid.nc

# because lacking of weighting file, I am going to use NCL to calcuate the hemispheric asymmetry.

###	# subsetting region and area average: SH
###	ncks -O -d lon,$lonS1,$lonE1 -d lat,$latS1,$latE1 ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid_SH.nc
###	ncwa -O -a lat,lon ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid_SH.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid_SH_avg.nc
###
###	# subsetting region and area average: NH
###	ncks -O -d lon,$lonS2,$lonE2 -d lat,$latS2,$latE2 ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid_NH.nc
###	ncwa -O -a lat,lon ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid_NH.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid_NH_avg.nc
###
###	# difference between SH and NH
###	ncdiff ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid_SH_avg.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid_NH_avg.nc ${dir[icase]}/${case[icase]}${append[icase]}_NetAtmHeat_regrid_SH-NH.nc
done
