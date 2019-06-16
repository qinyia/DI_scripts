
load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff_3D_2var(var_1[*][*][*]:numeric,lev[*]:numeric,lat_1[*]:numeric,lon_1[*]:numeric,\
		          var_2[*][*][*]:numeric,lat_2[*]:numeric,lon_2[*]:numeric)
begin
;horizontal interpolation
var_interp = lev_lat_lon_interpolate_2vars(var_1,lev,lat_1,lon_1, \
				           var_2,lat_2,lon_2)

var_interp_diff = var_interp(0,:,:,:)
var_interp_diff = var_interp(1,:,:,:) - var_interp(0,:,:,:)

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff_2D_2var(var_1[*][*]:numeric,lat_1[*]:numeric,lon_1[*]:numeric,\
		             var_2[*][*]:numeric,lat_2[*]:numeric,lon_2[*]:numeric)
begin
;horizontal interpolation
var_interp = lat_lon_interpolate_2vars(var_1,lat_1,lon_1, \
				       var_2,lat_2,lon_2)

var_interp_diff = var_interp(0,:,:)
var_interp_diff = var_interp(1,:,:) - var_interp(0,:,:)

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end

; %%%%%%%%%%%%%%%

function regrid_diff_2D_3var(var_obs[*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                          var_def[*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		          var_new[*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric,\
			  lat_spc[*]:numeric,lon_spc[*]:numeric)
begin
;horizontal interpolation
var_interp = lat_lon_interpolate_spcgrid(var_obs,lat_obs,lon_obs, \
				 var_def,lat_def,lon_def, \
				 var_new,lat_new,lon_new, \
				 lat_spc,lon_spc)

var_interp_diff = var_interp
var_interp_diff(0,:,:) = var_interp(1,:,:) - var_interp(0,:,:) ; def-obs
var_interp_diff(1,:,:) = var_interp(2,:,:) - var_interp(0,:,:) ; new-obs
var_interp_diff(2,:,:) = var_interp(2,:,:) - var_interp(1,:,:) ; new-obs

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
season					=(/"ANN"/)

var2D					=(/"U925","V925","TAUX","TAUY","SST","PRECT","CLDTOT","PSL"/)
compare_type 				= (/"cmip-obs"/)
case					= (/"def","new"/)

read_3D					= False
read_2D					= True

Dir					= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

if(compare_type(icomp) .eq."cmip-obs")then
	casename		 	= (/"","BC5_f19g16","BC5_f19g16_mac2"/)
	fdir 				= (/Dir+"obs/",Dir+"atm/",Dir+"atm/"/)
end if

if(compare_type(icomp) .eq."amip-obs")then
	casename		 	= (/"","FC5_f09f09_MG15_COSP","new-cldfrac16-1d-MG15-revise-2-COSP"/)
	fdir				= (/Dir+"obs/",Dir+"amip/",Dir+"amip/"/)
end if


;====================================
; deal with 2D variables
;====================================
if(read_2D)then

do is = 0,dimsizes(season)-1

fapend 					= (/"_"+season(is)+"_climo.nc", \
          	        		    "_"+season(is)+"_climo.nc", \
	  	        		    "_"+season(is)+"_climo.nc" /)

do ivar=0,dimsizes(var2D)-1

if(var2D(ivar).eq."SST")then
	casename(0) 			= "HadISST_PD"
	fdir 				= (/Dir+"obs/", \
	         			    Dir+"ocn/", \
		 			    Dir+"ocn/"/)
else
	fdir 				= (/Dir+"obs/", \
	         			    Dir+"atm/", \
		 			    Dir+"atm/"/)
end if

if(var2D(ivar).eq."PRECT")then
	casename(0) 			= "GPCP"
end if

if(var2D(ivar).eq."CLDLOW".or.var2D(ivar).eq."CLDMED".or.var2D(ivar).eq."CLDHGH".or.var2D(ivar).eq."CLDTOT")then
	casename(0) 			= "CLOUDSAT"
end if

if(var2D(ivar).eq."FSNS" .or. var2D(ivar).eq."FLNS" .or. \
   var2D(ivar).eq."LHFLX" .or. var2D(ivar).eq."SHFLX".or. \
   var2D(ivar).eq."TAUX".or.var2D(ivar).eq."TAUY".or. \
   var2D(ivar).eq."U925".or.var2D(ivar).eq."V925".or.\
   var2D(ivar).eq."PSL")then
	casename(0) 			= "JRA25"
end if


;;; for normal variables

result_obs				= get_var_2D(fdir(0),fapend(0),casename(0),var2D(ivar),"obs")
var_obs					= result_obs[0]
lat_obs					= result_obs[1]
lon_obs					= result_obs[2]

result_def				= get_var_2D(fdir(1),fapend(1),casename(1),var2D(ivar),"model")
var_def	 	       			= result_def[0]
lat_def	 	       			= result_def[1]
lon_def	 	       			= result_def[2]

result_new				= get_var_2D(fdir(2),fapend(2),casename(2),var2D(ivar),"model")
var_new					= result_new[0]
lat_new					= result_new[1]
lon_new					= result_new[2]

lat_spc					= fspan(-90,90,72)
lon_spc					= fspan(0,360,144)

result 					= regrid_diff_2D_3var(var_obs,lat_obs,lon_obs,\
	             				               var_def,lat_def,lon_def,\
		     				      	       var_new,lat_new,lon_new,\
							       lat_spc,lon_spc)
var_interp 				= result[0]
var_interp_diff 			= result[1]

delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])

printVarSummary(var_interp)

index					= dimsizes(var_interp)

if(ivar.eq.0)then
	vars_2D 			= new((/dimsizes(var2D),index(0),index(1),index(2)/),typeof(var_interp))
	vars_2D!0 			= "vars"
	vars_2D&vars 			= ispan(0,dimsizes(var2D)-1,1)
end if
vars_2D(ivar,:,:,:) 			= var_interp

delete(index)
delete(var_interp)
delete(var_interp_diff)


end do; do ivar=

printVarSummary(vars_2D)

;========================================================
; write into one NC file to check them
;========================================================
system("/bin/rm -f atm-var2D_"+compare_type(icomp)+"_"+season(is)+"_w_obs.nc")
ncdf = addfile("atm-var2D_"+compare_type(icomp)+"_"+season(is)+"_w_obs.nc" ,"c")

do ivar=0,dimsizes(var2D)-1

	ncdf->$var2D(ivar)$	 	= vars_2D(ivar,:,:,:)

end do ; do ivar=

;========================================================

end do; do is=

end if ; read_2D

end do; do icomp=
                               
end
