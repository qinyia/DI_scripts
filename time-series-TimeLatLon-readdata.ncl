
load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff_TimeLatLon_2var(var_1[*][*][*]:numeric,lat_1[*]:numeric,lon_1[*]:numeric,\
		                     var_2[*][*][*]:numeric,lat_2[*]:numeric,lon_2[*]:numeric)
begin
ntime		= dimsizes(var_1&time)

;horizontal interpolation
var_interp = time_lat_lon_interpolate_2vars(ntime,\
				            var_1,lat_1,lon_1, \
			                    var_2,lat_2,lon_2)

var_interp_diff = var_interp(0,:,:,:)
var_interp_diff = var_interp(1,:,:,:) - var_interp(0,:,:,:)

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin

;*********************************************
; read in data
;*********************************************

;season=(/"ANN","DJF","MAM","JJA","SON"/)
season			=(/"ANN"/)

;var_TimeLatLon		= (/"CLDLOW","ACRF","PRECT","TS","NetSfcHeat","PSL","FSNS","FLNS","LHFLX","SHFLX","NetAtmHeat","NetToaHeat","FSNT","FLNT","TAUX","TAUY"/)

var_TimeLatLon		=(/"TAUX","TAUY"/)

;var_TimeLatLon			= (/"NetSfcHeat","FSNS","TAUX","TAUY"/)

compare_type 		= (/"cmip-obs"/)
;compare_type 		= (/"amip-obs","cmip-obs"/)

compset			= (/"ocn"/)


do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

do icompset=0,dimsizes(compset)-1

if(compare_type(icomp) .eq."cmip-obs")then
Dir			= "/home/lyl/WORK1/lyl_cesm1_2_1/archive/"

	casename 	= (/"BC5_f19g16","BC5_f19g16_mac2"/)
	if(compset(icompset).eq."ocn")then
		fdir_append 	= (/"/ocn/hist/","/ocn/hist/"/)
	end if
	if(compset(icompset).eq."atm")then
		fdir_append	= (/"/atm/hist/","/atm/hist/"/)
	end if
end if

if(compare_type(icomp) .eq."amip-obs")then
Dir			= "/home/lyl/WORK4/cesm1_2_1/archive/"

	casename 	= (/"FC5_f19f19_MG15_cosp","FC5_f19f19_MG15_mac3_cosp"/)
	if(compset(icompset).eq."ocn")then
		print("")
		print("No available data for this compset:"+compset(icompset)+" in "+compare_type(icomp))
		print("")
	end if
	if(compset(icompset).eq."atm")then
		fdir_append	= (/"/atm/hist/","/atm/hist/"/)
	end if
end if

;====================================
; deal with 3D variables
;====================================

do ivar=0,dimsizes(var_TimeLatLon)-1

print(var_TimeLatLon(ivar))

do is = 0,dimsizes(season)-1

;;; for normal variables

if(compset(icompset).eq."atm")then
	result_def		= get_var_atm_TimeLatLon_addfiles(Dir,fdir_append(0),casename(0),var_TimeLatLon(ivar),"model")
	result_new		= get_var_atm_TimeLatLon_addfiles(Dir,fdir_append(1),casename(1),var_TimeLatLon(ivar),"model")
end if

if(compset(icompset).eq."ocn")then
	result_def		= get_var_ocn_TimeLatLon_addfiles(Dir,fdir_append(0),casename(0),var_TimeLatLon(ivar),"model")
	result_new		= get_var_ocn_TimeLatLon_addfiles(Dir,fdir_append(1),casename(1),var_TimeLatLon(ivar),"model")
end if

var_def	 	       	= result_def[0]
lat_def	 	       	= result_def[1]
lon_def	 	       	= result_def[2]
time_def		= result_def[3]

var_new			= result_new[0]
lat_new			= result_new[1]
lon_new			= result_new[2]
time_new		= result_new[3]

delete(result_def)
delete(result_new)



result 			= regrid_diff_TimeLatLon_2var(var_def,lat_def,lon_def,\
		             		              var_new,lat_new,lon_new)
var_interp 		= result[0]
var_interp_diff 	= result[1]


delete([/var_def,time_def,lat_def,lon_def,var_new,lat_new,lon_new,time_new/])

index			= dimsizes(var_interp)

if(ivar.eq.0)then
	vars_TimeLatLon 	= new((/dimsizes(var_TimeLatLon),index(0),index(1),index(2),index(3)/),typeof(var_interp))
	vars_TimeLatLon!0 	= "vars"
	vars_TimeLatLon&vars 	= ispan(0,dimsizes(var_TimeLatLon)-1,1)
end if
vars_TimeLatLon(ivar,:,:,:,:) 	= var_interp

delete(index)
delete(var_interp)
delete(var_interp_diff)


end do; do is=

end do; do ivar=

;========================================================
; write into one NC file to check them
;========================================================
system("/bin/rm -f vars_TimeLatLon_"+season+"_"+compset(icompset)+"_"+compare_type(icomp)+".nc")
ncdf = addfile("vars_TimeLatLon_"+season+"_"+compset(icompset)+"_"+compare_type(icomp)+".nc" ,"c")

;filedimdef(ncdf,"time",-1,True)

do ivar=0,dimsizes(var_TimeLatLon)-1
	ncdf->$var_TimeLatLon(ivar)$		= vars_TimeLatLon(ivar,:,:,:,:)
end do


end do ; do icompset=

end do; do icomp=
                               
end
