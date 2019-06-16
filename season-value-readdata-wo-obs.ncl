
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

function regrid_diff_time(time[*]:numeric,\
                          var_def[*][*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		          var_new[*][*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

;horizontal interpolation
var_interp = time_lat_lon_interpolate_2vars(dimsizes(time), \
					    var_def,lat_def,lon_def, \
					    var_new,lat_new,lon_new)


var_interp_diff = var_interp(0,:,:,:)
var_interp_diff = var_interp(1,:,:,:) - var_interp(0,:,:,:)

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin


;var_SeaLatLon=(/"PRECT","CLDTOT","CLDHGH","CLDMED","CLDLOW","FSNS","NetSfcHeat"/)
var_SeaLatLon=(/"PRECT","CLDTOT","CLDHGH","CLDMED","CLDLOW","FSNS","NetSfcHeat","TS","PSL","TAUX","TAUY"/)

;compare_type = (/"amip-obs","cmip-obs","cmip-amip"/)
;compare_type = (/"cmip-amip"/)
compare_type = (/"cmip-obs"/)

case=(/"def","new"/)

Dir="/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

do icase=0,dimsizes(case)-1
print(case(icase))

do ivar=0,dimsizes(var_SeaLatLon)-1


	if(compare_type(icomp) .eq."amip-obs")then
		casename = (/"FC5_f09f09_MG15_COSP","new-cldfrac16-1d-MG15-revise-2-COSP"/)
		fdir = (/Dir+"amip/",Dir+"amip/"/)

	end if
	if(compare_type(icomp) .eq."cmip-obs")then
		casename = (/"BC5_f19g16","BC5_f19g16_mac2"/)
		fdir = (/Dir+"atm/",Dir+"atm/"/)

	end if
	if(compare_type(icomp) .eq."cmip-amip")then
		fdir = (/Dir+"amip/",Dir+"atm/"/)
		if(case(icase).eq."def")then
			casename = (/"FC5_f09f09_MG15_COSP","BC5_f19g16"/)
		end if
		if(case(icase).eq."new")then
			casename = (/"new-cldfrac16-1d-MG15-revise-2-COSP","BC5_f19g16_mac2"/)
		end if
	end if



fapend = (/".cam.h0.MONTHS.nc", \
	   ".cam.h0.MONTHS.nc" /)


if(var_SeaLatLon(ivar).eq."SST")then
	fdir = (/Dir+"ocn/", \
		 Dir+"ocn/"/)
fapend = (/".pop.h.MONTHS.nc", \
	   ".pop.h.MONTHS.nc" /)

end if

result_def			= get_var_time(fdir(0),fapend(0),casename(0),var_SeaLatLon(ivar),"model")
var_def				= result_def[0]
lat_def				= result_def[1]
lon_def				= result_def[2]
time				= result_def[3]

result_new			= get_var_time(fdir(1),fapend(1),casename(1),var_SeaLatLon(ivar),"model")
var_new				= result_new[0]
lat_new				= result_new[1]
lon_new				= result_new[2]

printVarSummary(var_def)


result 				= regrid_diff_time(time, \
			 			   var_def,lat_def,lon_def,\
			 			   var_new,lat_new,lon_new)

var_interp 			= result[0]
var_interp_diff 		= result[1]

delete([/time,var_def,lat_def,lon_def,var_new,lat_new,lon_new,result/])


index					= dimsizes(var_interp)
if(ivar.eq.0)then
	vars_SeaLatLon			= new((/dimsizes(var_SeaLatLon),index(0),index(1),index(2),index(3)/),float)
	vars_SeaLatLon!0		= "vars"
	vars_SeaLatLon&vars		= ispan(0,dimsizes(var_SeaLatLon)-1,1)
end if

vars_SeaLatLon(ivar,:,:,:,:)		= var_interp


end do ; do ivar=

printVarSummary(vars_SeaLatLon)

;========================================================
; write into one NC file to check them
;========================================================
system("/bin/rm -f vars_SeaLatLon_"+compare_type(icomp)+"-wo-obs.nc")
ncdf = addfile("vars_SeaLatLon_"+compare_type(icomp)+"-wo-obs.nc","c")

do ivar=0,dimsizes(var_SeaLatLon)-1
	ncdf->$var_SeaLatLon(ivar)$		= vars_SeaLatLon(ivar,:,:,:,:)
end do

end do; do icase=

end do; do icomp=
 
end
