load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

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

function regrid_diff_3D_3var(var_obs[*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                          var_def[*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		          var_new[*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin
;horizontal interpolation
var_interp = lat_lon_interpolate(var_obs,lat_obs,lon_obs, \
				 var_def,lat_def,lon_def, \
				 var_new,lat_new,lon_new)

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

season			=(/"ANN"/)
var3D			=(/"TEMP","UVEL","VVEL","WVEL"/)
var2D			=(/"HMXL","SHF"/)
compare_type 		= (/"cmip-obs"/)
case			= (/"def","new"/)

Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

if(compare_type(icomp) .eq."cmip-obs")then
	casename 	= (/"","BC5_f19g16","BC5_f19g16_mac2"/)
	fdir 		= (/Dir+"obs/",Dir+"ocn/",Dir+"ocn/"/)
end if

;====================================
; deal with 3D variables
;====================================

do ivar=0,dimsizes(var3D)-1

do is = 0,dimsizes(season)-1

fapend 			= (/"_"+season(is)+"_climo.nc", \
          	            "_"+season(is)+"_climo.nc", \
	  	            "_"+season(is)+"_climo.nc" /)

;;; for normal variables

;result_obs		= get_var_ocn_3D_depth(fdir(0),fapend(0),casename(0),var3D(ivar),"obs")
;var_obs		= result_obs[0]
;lat_obs		= result_obs[1]
;lon_obs		= result_obs[2]

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


;result = regrid_diff(var_obs,lat_obs,lon_obs,\
;	             var_def,lat_def,lon_def,\
;		     var_new,lat_new,lon_new)
;var_interp = result[0]
;var_interp_diff = result[1]

result 			= regrid_diff_3D_2var(var_def,depth_def,lat_def,lon_def,\
		             		      var_new,lat_new,lon_new)
var_interp 		= result[0]
var_interp_diff 	= result[1]


;delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/var_def,depth_def,lat_def,lon_def,var_new,lat_new,lon_new,depth_new/])

index			= dimsizes(var_interp)

if(ivar.eq.0)then
	vars_3D 	= new((/dimsizes(var3D),index(0),index(1),index(2),index(3)/),typeof(var_interp))
	vars_3D!0 	= "vars"
	vars_3D&vars 	= var3D
end if
vars_3D(ivar,:,:,:,:) 	= var_interp

delete(index)
delete(var_interp)
delete(var_interp_diff)


end do; do is=

end do; do ivar=
printVarSummary(vars_3D)


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
	vars_2D&vars 	= var2D
end if
vars_2D(ivar,:,:,:) 	= var_interp

delete(index)
delete(var_interp)
delete(var_interp_diff)

end do; do is=

end do; do ivar=

printVarSummary(vars_2D)

;==============================================
; get the depth index of mixed layer depth
; get the averaged values over the mixed layer.
; get the entrainment velocity.
;==============================================

HMXL 			= 	vars_2D(0,:,:,:) ;unit is m.
depth 			= 	vars_3D&lev ; unit is m.

printVarSummary(HMXL)
printMinMax(HMXL,False)
printMinMax(depth,False)

nlat			= 	dimsizes(vars_2D&lat)
nlon			= 	dimsizes(vars_2D&lon)
ncase			=	dimsizes(vars_2D&case)

; for T, U, V, need to get the averaged values from sea surface to mixed layer depth.
; for W, should get the entrainment velocity, which is the vertical velocity at the mixed layer base.

index_mld 		= 	new((/ncase,nlat,nlon/),integer)
index_mld@longname	= 	"index of mixed layer depth in Depth dimension"
index_mld_tmp 		= 	index_mld

vars_3D_mldmean		=	new((/dimsizes(var3D),ncase,nlat,nlon/),float)

; entrainment velocity, the vertical velocity at the mixed layer base
w_mld			=	new((/ncase,nlat,nlon/),float)

; the temperature difference between the mixed layer and just below the mixed layer
; qinyi 2017-12-12 12:30:18 corrected
; the temperature difference should be the mixed layer mean temperature and the sublayer temperature

deltaT_mld_old		= 	new((/ncase,nlat,nlon/),float)
deltaT_mld		= 	new((/ncase,nlat,nlon/),float)

; for T, U, V, need to get the averaged values from sea surface to mixed layer depth.
; for W, should get the entrainment velocity, which is the vertical velocity at the mixed layer base.
do icase= 0,ncase-1
do ilat = 0,nlat-1
do ilon = 0,nlon-1

index_mld_tmp(icase,ilat,ilon) 		= 	closest_val(HMXL(icase,ilat,ilon),depth)	
if(.not.ismissing(index_mld_tmp(icase,ilat,ilon)).and.index_mld_tmp(icase,ilat,ilon).gt.1)then
	if(depth(index_mld_tmp(icase,ilat,ilon)).gt.HMXL(icase,ilat,ilon))then
		index_mld(icase,ilat,ilon)=index_mld_tmp(icase,ilat,ilon)-1
	else
		index_mld(icase,ilat,ilon)=index_mld_tmp(icase,ilat,ilon)
	end if
end if

if(.not.ismissing(index_mld(icase,ilat,ilon)))then

vars_3D_mldmean(:,icase,ilat,ilon) 	= dim_avg_n_Wrap(vars_3D(:,icase,:index_mld(icase,ilat,ilon),ilat,ilon),1)
w_mld(icase,ilat,ilon)			= (/vars_3D(3,icase,index_mld(icase,ilat,ilon),ilat,ilon)/)
deltaT_mld_old(icase,ilat,ilon)		= (/vars_3D(0,icase,index_mld(icase,ilat,ilon),ilat,ilon)-vars_3D(0,icase,index_mld(icase,ilat,ilon)+2,ilat,ilon)/)
deltaT_mld(icase,ilat,ilon)		= (/vars_3D_mldmean(0,icase,ilat,ilon)-vars_3D(0,icase,index_mld(icase,ilat,ilon)+2,ilat,ilon)/)

else

vars_3D_mldmean(:,icase,ilat,ilon) 	=	vars_3D@_FillValue
w_mld(icase,ilat,ilon)			= 	vars_3D@_FillValue
deltaT_mld_old(icase,ilat,ilon)		=	vars_3D@_FillValue
deltaT_mld(icase,ilat,ilon)		=	vars_3D@_FillValue

end if

end do ; do ilon=
end do ; do ilat=
end do ; do icase=

vars_3D_mldmean!0			= "vars"
vars_3D_mldmean&vars			= vars_3D&vars
vars_3D_mldmean!1			= "case"
vars_3D_mldmean&case			= vars_3D&case
vars_3D_mldmean!2			= "lat"
vars_3D_mldmean&lat			= vars_3D&lat
vars_3D_mldmean!3			= "lon"
vars_3D_mldmean&lon			= vars_3D&lon

copy_VarCoords_1(vars_3D_mldmean(case|:,lat|:,lon|:,vars|:),w_mld)
copy_VarCoords_1(vars_3D_mldmean(case|:,lat|:,lon|:,vars|:),deltaT_mld_old)
copy_VarCoords_1(vars_3D_mldmean(case|:,lat|:,lon|:,vars|:),deltaT_mld)
copy_VarCoords_1(vars_3D_mldmean(case|:,lat|:,lon|:,vars|:),index_mld)

printVarSummary(index_mld)
printVarSummary(vars_3D_mldmean)
printVarSummary(w_mld)
printMinMax(w_mld,False)
printVarSummary(deltaT_mld)
printMinMax(deltaT_mld_old,True)
printMinMax(deltaT_mld,True)

;========================================================
; calculate entrainment term
;========================================================

entrain			= HMXL
entrain 		= w_mld*deltaT_mld/HMXL
entrain			= entrain*(-1.0)

entrain 		= entrain*30*24*3600.	;change from K/s to K/month

entrain@longname	= "entrainment term"	
entrain@units		= "degC/month"

printVarSummary(entrain)
printMinMax(entrain(:,{-20:0},{180:240}),False)

;========================================================
; calculate horizontal advection: zonal and meridional
;========================================================
; get latitudinal and longitudinal gradient
; here, leave one small problem: how to directly use "TEMP" as the index?
; temperature is the first variable in var3D, so its index is 0.
Temp_mldmean 		= 	vars_3D_mldmean(0,:,:,:)
;Temp_mldmean 		= 	vars_3D(0,:,3,:,:)


;---------------------------------------------------------------
; use center_finite_diff_n to test the "wave-like" structure
; use forward finite diff to test the structure
; still not work....
;---------------------------------------------------------------

;lon			= 	vars_3D_mldmean&lon
;dlon			=	(lon(2)-lon(1))*0.0174533
;
;lat			= 	vars_3D_mldmean&lat
;dlat			=	(lat(2)-lat(1))*0.0174533
;
;T_grad_dy		= Temp_mldmean
;T_grad_dx		= Temp_mldmean
;
;do nl=0,nlat-1
;do ll=0,nlon-2
;	dX 			= 	6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
;
;;	T_grad_dx(:,nl:nl,:) 	= 	center_finite_diff_n (Temp_mldmean(:,nl:nl,:), dX , True,0,2)
;	T_grad_dx(:,nl,ll) 	= (Temp_mldmean(:,nl,ll+1)-Temp_mldmean(:,nl,ll))/dX
;end do
;T_grad_dx(:,nl,nlon-1)		= (Temp_mldmean(:,nl,nlon-1)-Temp_mldmean(:,nl,nlon-2))/dX
;end do
;
;do nl=0,nlon-1
;	dY 			= 	6378388.*dlat
;	T_grad_dy(:,:,nl:nl) 	= 	center_finite_diff_n (Temp_mldmean(:,:,nl:nl), dY , True,0,1)
;end do

;---------------------------------------------------------------


gradLatLon  		= 	grad_latlon_cfd (Temp_mldmean, vars_3D&lat, vars_3D&lon, True, False)
                	               
T_grad_dy  		= 	gradLatLon[0] 
T_grad_dx  		= 	gradLatLon[1]

delete(gradLatLon)

T_grad_dx@long_name 	= 	"zonal gradient"
T_grad_dy@long_name 	= 	"meridional gradient"
T_grad_dx@units     	= 	"degC/m"
T_grad_dy@units     	= 	"degC/m"

; zonal advection: u_dTdx
;u_mldmean 		= 	vars_3D(1,:,3,:,:)
u_mldmean		= 	vars_3D_mldmean(1,:,:,:)

u_dTdx 			= 	T_grad_dx
u_dTdx 			= 	u_mldmean*T_grad_dx
u_dTdx			=	u_dTdx*(-1.)

;change unit from K/s to K/month
u_dTdx 			= 	u_dTdx*30*24*3600

u_dTdx@long_name	= 	"advection term: u_dTdx"
u_dTdx@units 		= 	"degC/month"

printVarSummary(u_dTdx)
printMinMax(u_dTdx,False)

; meridional advection: v_dTdy
;v_mldmean			= 	vars_3D(2,:,3,:,:)
v_mldmean		= 	vars_3D_mldmean(2,:,:,:)

v_dTdy 			= 	T_grad_dy
v_dTdy 			= 	v_mldmean*T_grad_dy
v_dTdy			= 	v_dTdy*(-1.)

;change units
v_dTdy 			= 	v_dTdy*30*24*3600

v_dTdy@longname 	= 	"advection term: v_dTdy"
v_dTdy@units 		= 	"degC/month"

printVarSummary(v_dTdy)
printMinMax(v_dTdy,False)


;========================================================
; calculate the net energy flux term
;========================================================

rho_sw 			= 1030.		; kg/m3
cp_sw			= 4.0e3		; J/kg/K

SHF			= vars_2D(1,:,:,:) ; total surface heat flux: (case,lat,lon)
printVarSummary(SHF)
printMinMax(SHF,False)

netFlux 		= SHF/rho_sw/cp_sw/HMXL
copy_VarCoords(SHF,netFlux)

netFlux			= netFlux*30*24*3600.
netFlux@longname 	= "net energy flux term"
netFlux@units		= "degC/month"

printVarSummary(netFlux)
printMinMax(netFlux(:,{-20:0},{180:240}),False)

;========================================================
; calculate the sum of four terms:
; net energy flux: netFlux
; zonal advection: u_dTdx
; meridional advection: v_dTdy
; entrainment: entrain
;========================================================

sumall 			= netFlux + u_dTdx + v_dTdy + entrain
copy_VarCoords(netFlux,sumall)
residual		= (-1.)*sumall
copy_VarCoords(netFlux,residual)
sumall@longname 	= "sum of four terms"
sumall@units 		= "degC/month"

residual@longname 	= "residual"
residual@units 		= "degC/month"

printVarSummary(sumall)


;========================================================
; change units from degC/month to W/m2
;========================================================

factor			= 30.*24.*3600.
residual_Wm2		= residual*cp_sw*rho_sw*HMXL/factor
sumall_Wm2		= sumall*cp_sw*rho_sw*HMXL/factor
netFlux_Wm2		= netFlux*cp_sw*rho_sw*HMXL/factor
u_dTdx_Wm2		= u_dTdx*cp_sw*rho_sw*HMXL/factor
v_dTdy_Wm2		= v_dTdy*cp_sw*rho_sw*HMXL/factor
entrain_Wm2		= entrain*cp_sw*rho_sw*HMXL/factor

copy_VarCoords(residual,residual_Wm2)
copy_VarCoords(sumall,sumall_Wm2)
copy_VarCoords(netFlux,netFlux_Wm2)
copy_VarCoords(u_dTdx,u_dTdx_Wm2)
copy_VarCoords(v_dTdy,v_dTdy_Wm2)
copy_VarCoords(entrain,entrain_Wm2)


;========================================================
; put all terms into one array to faciliate the plots
;========================================================
output_terms		= (/"residual","netFlux","u_dTdx","v_dTdy","entrain"/)
terms			= ispan(0,dimsizes(output_terms)-1,1)
finalvar		= new((/dimsizes(output_terms),ncase,nlat,nlon/),float)
finalvar!0		= "terms"
finalvar&terms		= terms

finalvar(0,:,:,:)	= netFlux
finalvar		= (/residual,netFlux,u_dTdx,v_dTdy,entrain/)


finalvardiff 		= new((/dimsizes(output_terms),nlat,nlon/),float)
finalvardiff!0		= "terms"
finalvardiff&terms	= terms
finalvardiff!1		= "lat"
finalvardiff&lat	= netFlux&lat
finalvardiff!2		= "lon"
finalvardiff&lon	= netFlux&lon

finalvardiff		= (/residual(1,:,:)-residual(0,:,:),\
			    netFlux(1,:,:)-netFlux(0,:,:),\
			    u_dTdx(1,:,:)-u_dTdx(0,:,:),\
			    v_dTdy(1,:,:)-v_dTdy(0,:,:),\
			    entrain(1,:,:)-entrain(0,:,:)/)


finalvar_Wm2		= finalvar
finalvar_Wm2		= (/residual_Wm2,netFlux_Wm2,u_dTdx_Wm2,v_dTdy_Wm2,entrain_Wm2/)

finalvardiff_Wm2 	= finalvardiff
finalvardiff_Wm2	= (/residual_Wm2(1,:,:)-residual_Wm2(0,:,:),\
			    netFlux_Wm2(1,:,:)-netFlux_Wm2(0,:,:),\
			    u_dTdx_Wm2(1,:,:)-u_dTdx_Wm2(0,:,:),\
			    v_dTdy_Wm2(1,:,:)-v_dTdy_Wm2(0,:,:),\
			    entrain_Wm2(1,:,:)-entrain_Wm2(0,:,:)/)

printVarSummary(finalvar_Wm2)
printMinMax(finalvar_Wm2,False)

printVarSummary(finalvardiff_Wm2)
printMinMax(finalvardiff_Wm2,False)

;====================================================================
; output the difference of some basic variables between two schemes
;====================================================================

T_grad_dx_diff		= T_grad_dx(1,:,:)-T_grad_dx(0,:,:)
copy_VarCoords_1(T_grad_dx(lat|:,lon|:,case|:),T_grad_dx)

T_grad_dy_diff		= T_grad_dy(1,:,:)-T_grad_dy(0,:,:)
copy_VarCoords_1(T_grad_dy(lat|:,lon|:,case|:),T_grad_dy)

printVarSummary(T_grad_dy_diff)
printMinMax(T_grad_dy_diff,True)

Temp_mldmean_diff	= Temp_mldmean(1,:,:)-Temp_mldmean(0,:,:)
copy_VarCoords_1(Temp_mldmean(lat|:,lon|:,case|:),Temp_mldmean_diff)

printVarSummary(Temp_mldmean_diff)
printMinMax(Temp_mldmean_diff,True)

u_mldmean_diff		= u_mldmean(1,:,:)-u_mldmean(0,:,:)
copy_VarCoords_1(u_mldmean(lat|:,lon|:,case|:),u_mldmean_diff)

v_mldmean_diff		= v_mldmean(1,:,:)-v_mldmean(0,:,:)
copy_VarCoords_1(v_mldmean(lat|:,lon|:,case|:),v_mldmean_diff)

printVarSummary(v_mldmean_diff)
printMinMax(v_mldmean_diff,True)

w_mld_diff		= w_mld(1,:,:)-w_mld(0,:,:)
copy_VarCoords_1(w_mld(lat|:,lon|:,case|:),w_mld_diff)

printVarSummary(w_mld_diff)
printMinMax(w_mld_diff,True)

deltaT_mld_diff		= deltaT_mld(1,:,:)-deltaT_mld(0,:,:)
copy_VarCoords_1(deltaT_mld(lat|:,lon|:,case|:),deltaT_mld_diff)

printVarSummary(deltaT_mld_diff)
printMinMax(deltaT_mld_diff,True)

HMXL_diff		= HMXL(1,:,:)-HMXL(0,:,:)
copy_VarCoords_1(HMXL(lat|:,lon|:,case|:),HMXL_diff)

printVarSummary(HMXL_diff)
printMinMax(HMXL_diff,True)

residual_Wm2_diff	= residual_Wm2(1,:,:)-residual_Wm2(0,:,:)
copy_VarCoords_1(residual_Wm2(lat|:,lon|:,case|:),residual_Wm2_diff)

sumall_Wm2_diff		= sumall_Wm2(1,:,:)-sumall_Wm2(0,:,:)
copy_VarCoords_1(sumall_Wm2(lat|:,lon|:,case|:),sumall_Wm2_diff)

netFlux_Wm2_diff	= netFlux_Wm2(1,:,:)-netFlux_Wm2(0,:,:)
copy_VarCoords_1(netFlux_Wm2(lat|:,lon|:,case|:),netFlux_Wm2_diff)

u_dTdx_Wm2_diff		= u_dTdx_Wm2(1,:,:)-u_dTdx_Wm2(0,:,:)
copy_VarCoords_1(u_dTdx_Wm2(lat|:,lon|:,case|:),u_dTdx_Wm2_diff)

v_dTdy_Wm2_diff		= v_dTdy_Wm2(1,:,:)-v_dTdy_Wm2(0,:,:)
copy_VarCoords_1(v_dTdy_Wm2(lat|:,lon|:,case|:),v_dTdy_Wm2_diff)

entrain_Wm2_diff	= entrain_Wm2(1,:,:)-entrain_Wm2(0,:,:)
copy_VarCoords_1(entrain_Wm2(lat|:,lon|:,case|:),entrain_Wm2_diff)

;========================================================
; write into one NC file to check them
;========================================================
system("/bin/rm -f terms"+"_"+season+".nc")
ncdf = addfile("terms"+"_"+season+".nc" ,"c")

;filedimdef(ncdf,"time",-1,True)

ncdf->residual_Wm2 	= residual_Wm2
ncdf->sumall_Wm2 	= sumall_Wm2
ncdf->netFlux_Wm2	= netFlux_Wm2
ncdf->u_dTdx_Wm2 	= u_dTdx_Wm2
ncdf->v_dTdy_Wm2 	= v_dTdy_Wm2
ncdf->entrain_Wm2	= entrain_Wm2

ncdf->residual 	= residual
ncdf->sumall 	= sumall
ncdf->netFlux	= netFlux
ncdf->u_dTdx 	= u_dTdx
ncdf->v_dTdy 	= v_dTdy
ncdf->entrain	= entrain

ncdf->residual_Wm2_diff = residual_Wm2_diff
ncdf->sumall_Wm2_diff 	= sumall_Wm2_diff
ncdf->netFlux_Wm2_diff	= netFlux_Wm2_diff
ncdf->u_dTdx_Wm2_diff	= u_dTdx_Wm2_diff
ncdf->v_dTdy_Wm2_diff	= v_dTdy_Wm2_diff
ncdf->entrain_Wm2_diff	= entrain_Wm2_diff

ncdf->finalvar_Wm2	= finalvar_Wm2
ncdf->finalvardiff_Wm2	= finalvardiff_Wm2

ncdf->T_grad_dx		= T_grad_dx
ncdf->T_grad_dy		= T_grad_dy
ncdf->Temp_mldmean	= Temp_mldmean
ncdf->u_mldmean		= u_mldmean
ncdf->v_mldmean		= v_mldmean
ncdf->w_mld		= w_mld
ncdf->deltaT_mld	= deltaT_mld
ncdf->deltaT_mld_old	= deltaT_mld_old
ncdf->HMXL		= HMXL

ncdf->T_grad_dx_diff	= T_grad_dx_diff
ncdf->T_grad_dy_diff	= T_grad_dy_diff
ncdf->Temp_mldmean_diff	= Temp_mldmean_diff
ncdf->u_mldmean_diff	= u_mldmean_diff
ncdf->v_mldmean_diff	= v_mldmean_diff
ncdf->w_mld_diff	= w_mld_diff
ncdf->deltaT_mld_diff	= deltaT_mld_diff
ncdf->HMXL_diff		= HMXL_diff

ncdf->index_mld		= index_mld

end do; do icomp=
                               
end
