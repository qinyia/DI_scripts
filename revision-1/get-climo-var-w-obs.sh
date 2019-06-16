
dir1=/home/lyl/WORK3/qinyi/double-ITCZ/data/obs/
dir2=/home/lyl/WORK3/qinyi/double-ITCZ/data/atm/

case=(obs BC5_f19g16 BC5_f19g16_mac2)
ncase=${#case[@]}

#obscase=(CERES-EBAF CERES-EBAF ISCCP ISCCP ISCCP ISCCP CERES-EBAF CERES-EBAF)
obscase=(CERES-EBAF CERES-EBAF CLOUDSAT CLOUDSAT CLOUDSAT CLOUDSAT CERES-EBAF CERES-EBAF)
obsvar=(FSNTOA FLUT CLDHGH CLDMED CLDLOW CLDTOT SWCF LWCF)
nobsvar=${#obsvar[@]}

modvar=(FSNT,FLNT,CLDHGH,CLDMED,CLDLOW,CLDTOT,SWCF,LWCF)

append=(_ANN_climo _ANN_climo _ANN_climo)
dir=($dir1 $dir2 $dir2)
var=($var1 $var2 $var2)



workdir=/home/lyl/WORK3/qinyi/double-ITCZ/
cd $workdir


for icase in `seq 0 $[$ncase-1]`
do

rm ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc

	if [ "$icase" = "0" ];then
		# extract obs cloud variables from dataset
		for ivar in `seq 0 $[$nobsvar-1]`
		do
			# extract each var from each obs data
			if [ "${obscase[ivar]}" == "CLOUDSAT" ];then 
				#add time dimension
                		ncecat -O -u time ${dir[icase]}/${obscase[ivar]}${append[icase]}.nc ${dir[icase]}/${obscase[ivar]}${append[icase]}_w_timedim.nc
                		#add time variable
                		ncap2 -O -s 'time=array(1,1,$time)' ${dir[icase]}/${obscase[ivar]}${append[icase]}_w_timedim.nc ${dir[icase]}/${obscase[ivar]}${append[icase]}_w_timedim.nc
				ncks -O -v ${obsvar[ivar]} ${dir[icase]}/${obscase[ivar]}${append[icase]}_w_timedim.nc ${dir[icase]}/${obsvar[ivar]}.nc
			else
				ncks -O -v ${obsvar[ivar]} ${dir[icase]}/${obscase[ivar]}${append[icase]}.nc ${dir[icase]}/${obsvar[ivar]}.nc
			fi
			# regrid data
			cdo remapbil,regrid.txt ${dir[icase]}/${obsvar[ivar]}.nc ${dir[icase]}/${obsvar[ivar]}_regrid.nc
			rm ${dir[icase]}/${obsvar[ivar]}.nc
			# append all obs variables in one file
			ncks -A ${dir[icase]}/${obsvar[ivar]}_regrid.nc ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc
			rm ${dir[icase]}/${obsvar[ivar]}_regrid.nc
		done
			# rename variable
			ncrename -O -v FSNTOA,FSNT ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc
			ncrename -O -v FLUT,FLNT ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc
			# calculate net TOA heat flux
			ncap2 -O -s 'NetTOAHeat=(FSNT-FLNT)' ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc
			exit

	else
		# extract model cloud variables from data
		ncks -O -v ${modvar} ${dir[icase]}/${case[icase]}${append[icase]}.nc ${dir[icase]}/${case[icase]}_glbvar_for_glbmean.nc
		# regrid
		cdo remapbil,regrid.txt ${dir[icase]}/${case[icase]}_glbvar_for_glbmean.nc ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc
			# calculate net TOA heat flux
			ncap2 -O -s 'NetTOAHeat=(FSNT-FLNT)' ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc
		# change cloud fraction units to percent
		ncap2 -O -s 'CLDLOW=CLDLOW*100' ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc  ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc
		ncap2 -O -s 'CLDMED=CLDMED*100' ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc  ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc
		ncap2 -O -s 'CLDHGH=CLDHGH*100' ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc  ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc
		ncap2 -O -s 'CLDTOT=CLDTOT*100' ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc  ${dir[icase]}/${case[icase]}_glbvar_for_glbmean_regrid.nc

	fi

done
