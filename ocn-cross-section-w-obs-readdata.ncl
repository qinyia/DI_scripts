
load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

; %%%%%%%%%%%%%%%

function regrid_diff_3D_3var(var_obs[*][*][*]:numeric,depth_obs[*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                             var_def[*][*][*]:numeric,depth_def[*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		             var_new[*][*][*]:numeric,depth_new[*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

; vertical interpolation: Interpolate to observation depths
var_def_tmp		= int2p_n_Wrap(depth_def,var_def,depth_obs,1,0)
var_new_tmp		= int2p_n_Wrap(depth_new,var_new,depth_obs,1,0)
printVarSummary(var_def_tmp)

;horizontal interpolation
var_interp = lev_lat_lon_interpolate(var_obs,depth_obs,lat_obs,lon_obs, \
				     var_def_tmp,lat_def,lon_def, \
				     var_new_tmp,lat_new,lon_new)

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

;*********************************************
; read in data
;*********************************************

;season=(/"ANN","DJF","MAM","JJA","SON"/)
season			=(/"ANN"/)
var3D			=(/"TEMP","UVEL","VVEL"/)

var2D			=(/"TAUX","TAUY"/)
compare_type 		= (/"cmip-obs"/)
case			= (/"def","new"/)

Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

if(compare_type(icomp) .eq."cmip-obs")then
	casename 	= (/"SODA","BC5_f19g16","BC5_f19g16_mac2"/)
	fdir 		= (/Dir+"obs/",Dir+"ocn/",Dir+"ocn/"/)
end if

;====================================
; deal with 3D variables
;====================================


do is = 0,dimsizes(season)-1

fapend 			= (/"_"+season(is)+"_climo.nc", \
          	            "_"+season(is)+"_climo.nc", \
	  	            "_"+season(is)+"_climo.nc" /)

do ivar=0,dimsizes(var3D)-1
;;; for normal variables

result_obs		= get_var_ocn_3D_depth(fdir(0),fapend(0),casename(0),var3D(ivar),"obs")
var_obs			= result_obs[0]
lat_obs			= result_obs[1]
lon_obs			= result_obs[2]
depth_obs		= result_obs[3]
printVarSummary(var_obs)
printMinMax(var_obs,True)

result_def		= get_var_ocn_3D_depth(fdir(1),fapend(1),casename(1),var3D(ivar),"model")
var_def	 	       	= result_def[0]
lat_def	 	       	= result_def[1]
lon_def	 	       	= result_def[2]
depth_def		= result_def[3]

result_new		= get_var_ocn_3D_depth(fdir(2),fapend(2),casename(2),var3D(ivar),"model")
var_new			= result_new[0]
lat_new			= result_new[1]
lon_new			= result_new[2]
depth_new		= result_new[3]


result 			= regrid_diff_3D_3var(var_obs,depth_obs,lat_obs,lon_obs,\
					      var_def,depth_def,lat_def,lon_def,\
		             		      var_new,depth_new,lat_new,lon_new)
var_interp 		= result[0]
var_interp_diff 	= result[1]
printVarSummary(var_interp)

delete([/var_obs,depth_obs,lat_obs,lon_obs,var_def,depth_def,lat_def,lon_def,var_new,depth_new,lat_new,lon_new/])

index			= dimsizes(var_interp)

if(ivar.eq.0)then
	vars_3D 	= new((/dimsizes(var3D),index(0),index(1),index(2),index(3)/),typeof(var_interp))
	vars_3D!0 	= "vars"
	vars_3D&vars 	= ispan(0,dimsizes(var3D)-1,1)
end if
vars_3D(ivar,:,:,:,:) 	= var_interp

delete(index)
delete(var_interp)
delete(var_interp_diff)


end do; do ivar=

printVarSummary(vars_3D)

;========================================================
; write into one NC file to check them
;========================================================

system("/bin/rm -f ocn-cross-section-w-obs-vars_3D_"+season(is)+".nc")
ncdf = addfile("ocn-cross-section-w-obs-vars_3D_"+season(is)+".nc" ,"c")

do ivar=0,dimsizes(var3D)-1
ncdf->$var3D(ivar)$	 	= vars_3D(ivar,:,:,:,:)
end do

end do; do is=

exit



;====================================
; deal with 2D variables
;====================================

do ivar=0,dimsizes(var2D)-1


do is = 0,dimsizes(season)-1

fapend 			= (/"_"+season(is)+"_climo.nc", \
       			    "_"+season(is)+"_climo.nc", \
		            "_"+season(is)+"_climo.nc" /)

;;; for normal variables

;result_obs		= get_var_ocn_2D(fdir(0),fapend(0),casename(0),var2D(ivar),"obs")
;var_obs		= result_obs[0]
;lat_obs		= result_obs[1]
;lon_obs		= result_obs[2]

result_def		= get_var_ocn_2D(fdir(1),fapend(1),casename(1),var2D(ivar),"model")
var_def			= result_def[0]
lat_def			= result_def[1]
lon_def			= result_def[2]

result_new		= get_var_ocn_2D(fdir(2),fapend(2),casename(2),var2D(ivar),"model")
var_new			= result_new[0]
lat_new			= result_new[1]
lon_new			= result_new[2]


;result 		= regrid_diff(var_obs,lat_obs,lon_obs,\
;	            	 	      var_def,lat_def,lon_def,\
;		    	              var_new,lat_new,lon_new)
;var_interp 		= result[0]
;var_interp_diff 	= result[1]

result 			= regrid_diff_2D_2var(var_def,lat_def,lon_def,\
		                      var_new,lat_new,lon_new)
var_interp 		= result[0]
var_interp_diff 	= result[1]


;delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/var_def,lat_def,lon_def,var_new,lat_new,lon_new/])

index			= dimsizes(var_interp)

if(ivar.eq.0)then
	vars_2D 	= new((/dimsizes(var2D),index(0),index(1),index(2)/),typeof(var_interp))
	vars_2D!0 	= "vars"
	vars_2D&vars 	= ispan(0,dimsizes(var2D)-1,1)
end if
vars_2D(ivar,:,:,:) 	= var_interp

delete(index)
delete(var_interp)
delete(var_interp_diff)

end do; do is=

end do; do ivar=

printVarSummary(vars_2D)

;========================================================
; write into one NC file to check them
;========================================================
system("/bin/rm -f vars_2D_"+season+".nc")
ncdf = addfile("vars_2D_"+season+".nc" ,"c")

ncdf->vars_2D	 	= vars_2D


end do; do icomp=
                               
end
