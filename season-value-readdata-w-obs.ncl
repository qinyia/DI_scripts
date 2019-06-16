
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

function regrid_diff_time(time[*]:numeric,\
                          var_obs[*][*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                          var_def[*][*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		          var_new[*][*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric,\
			  lat_spc[*]:numeric,lon_spc[*]:numeric)
begin

;horizontal interpolation
var_interp = time_lat_lon_interpolate_3vars_specgrid(dimsizes(time), \
                                           	     var_obs,lat_obs,lon_obs, \
					   	     var_def,lat_def,lon_def, \
					   	     var_new,lat_new,lon_new, \
						     lat_spc,lon_spc)

var_interp_diff = var_interp
var_interp_diff(0,:,:,:) = var_interp(1,:,:,:) - var_interp(0,:,:,:) ; def-obs
var_interp_diff(1,:,:,:) = var_interp(2,:,:,:) - var_interp(0,:,:,:) ; new-obs
var_interp_diff(2,:,:,:) = var_interp(2,:,:,:) - var_interp(1,:,:,:) ; new-obs

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin


var_SeaLatLon=(/"SST","PRECT","CLDTOT","CLDHGH","CLDMED","CLDLOW","FSNS","NetSfcHeat"/)
;var_SeaLatLon=(/"NetSfcHeat"/)
;var_SeaLatLon=(/"FSNS","FLNS","LHFLX","SHFLX"/)
;var_SeaLatLon=(/"FLNS"/)
;var_SeaLatLon=(/"LHFLX"/)
;var_SeaLatLon=(/"CLDTOT"/)

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
		casename = (/"","FC5_f09f09_MG15_COSP","new-cldfrac16-1d-MG15-revise-2-COSP"/)
		fdir = (/Dir+"obs/",Dir+"amip/",Dir+"amip/"/)

	end if
	if(compare_type(icomp) .eq."cmip-obs")then
		casename = (/"","BC5_f19g16","BC5_f19g16_mac2"/)
		fdir = (/Dir+"obs/",Dir+"atm/",Dir+"atm/"/)

	end if
	if(compare_type(icomp) .eq."cmip-amip")then
		fdir = (/Dir+"obs/",Dir+"amip/",Dir+"atm/"/)
		if(case(icase).eq."def")then
			casename = (/"","FC5_f09f09_MG15_COSP","BC5_f19g16"/)
		end if
		if(case(icase).eq."new")then
			casename = (/"","new-cldfrac16-1d-MG15-revise-2-COSP","BC5_f19g16_mac2"/)
		end if
	end if



fapend = (/"_MONTHS.nc", \
           ".cam.h0.MONTHS.nc", \
	   ".cam.h0.MONTHS.nc" /)


if(var_SeaLatLon(ivar).eq."SST")then
	casename(0) = "HadISST_PD"
	fdir = (/Dir+"obs/", \
	         Dir+"ocn/", \
		 Dir+"ocn/"/)
fapend = (/"_MONTHS.nc", \
           ".pop.h.MONTHS.nc", \
	   ".pop.h.MONTHS.nc" /)

end if

if(var_SeaLatLon(ivar).eq."PRECT")then
	casename(0) = "GPCP"
end if

if(var_SeaLatLon(ivar).eq."CLDLOW".or.var_SeaLatLon(ivar).eq."CLDMED".or.var_SeaLatLon(ivar).eq."CLDHGH".or.var_SeaLatLon(ivar).eq."CLDTOT")then
	casename(0) = "CLOUDSAT"
end if

if(var_SeaLatLon(ivar).eq."FSNS" .or. var_SeaLatLon(ivar).eq."FLNS" .or. \
   var_SeaLatLon(ivar).eq."LHFLX" .or. var_SeaLatLon(ivar).eq."SHFLX" .or.\
   var_SeaLatLon(ivar).eq."NetSfcHeat")then
	casename(0) = "JRA25"
end if



result_obs			= get_var_time(fdir(0),fapend(0),casename(0),var_SeaLatLon(ivar),"obs")
var_obs				= result_obs[0]
lat_obs				= result_obs[1]
lon_obs				= result_obs[2]
time				= result_obs[3]

result_def			= get_var_time(fdir(1),fapend(1),casename(1),var_SeaLatLon(ivar),"model")
var_def				= result_def[0]
lat_def				= result_def[1]
lon_def				= result_def[2]

result_new			= get_var_time(fdir(2),fapend(2),casename(2),var_SeaLatLon(ivar),"model")
var_new				= result_new[0]
lat_new				= result_new[1]
lon_new				= result_new[2]

printVarSummary(var_obs)

lat_spc				= fspan(-90,90,72)
lon_spc				= fspan(0,360,144)

result 				= regrid_diff_time(time, \
			 			   var_obs,lat_obs,lon_obs,\
			 			   var_def,lat_def,lon_def,\
			 			   var_new,lat_new,lon_new,\
						   lat_spc,lon_spc)

var_interp 			= result[0]
var_interp_diff 		= result[1]

delete([/var_obs,lat_obs,lon_obs,time,var_def,lat_def,lon_def,var_new,lat_new,lon_new,result/])


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
system("/bin/rm -f vars_SeaLatLon_"+compare_type(icomp)+"-w-obs.nc")
ncdf = addfile("vars_SeaLatLon_"+compare_type(icomp)+"-w-obs.nc","c")

do ivar=0,dimsizes(var_SeaLatLon)-1
	ncdf->$var_SeaLatLon(ivar)$		= vars_SeaLatLon(ivar,:,:,:,:)
end do

end do; do icase=

end do; do icomp=
 
end
