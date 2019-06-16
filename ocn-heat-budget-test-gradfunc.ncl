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
season=(/"ANN"/)

var3D=(/"TEMP","UVEL","VVEL","WVEL"/)
var2D=(/"HMXL","SHF"/)


compare_type = (/"cmip-obs"/)

case=(/"def","new"/)

Dir="/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

do icase=0,dimsizes(case)-1

	if(compare_type(icomp) .eq."cmip-obs")then
		casename = (/"","BC5_f19g16","BC5_f19g16_mac2"/)
		fdir = (/Dir+"obs/",Dir+"ocn/",Dir+"ocn/"/)
	end if

;====================================
; deal with 3D variables
;====================================

do ivar=0,dimsizes(var3D)-1

do is = 0,dimsizes(season)-1

fapend = (/"_"+season(is)+"_climo.nc", \
           "_"+season(is)+"_climo.nc", \
	   "_"+season(is)+"_climo.nc" /)

;;; for normal variables

;result_obs	= get_var_ocn_3D_depth(fdir(0),fapend(0),casename(0),var3D(ivar),"obs")
;var_obs		= result_obs[0]
;lat_obs		= result_obs[1]
;lon_obs		= result_obs[2]

result_def	= get_var_ocn_3D_depth(fdir(1),fapend(1),casename(1),var3D(ivar),"model")
var_def		= result_def[0]
lat_def		= result_def[1]
lon_def		= result_def[2]
depth_def	= result_def[3]

result_new	= get_var_ocn_3D_depth(fdir(2),fapend(2),casename(2),var3D(ivar),"model")
var_new		= result_new[0]
lat_new		= result_new[1]
lon_new		= result_new[2]
depth_new	= result_new[3]


;result = regrid_diff(var_obs,lat_obs,lon_obs,\
;	             var_def,lat_def,lon_def,\
;		     var_new,lat_new,lon_new)
;var_interp = result[0]
;var_interp_diff = result[1]

result = regrid_diff_3D_2var(var_def,depth_def,lat_def,lon_def,\
		             var_new,lat_new,lon_new)
var_interp = result[0]
var_interp_diff = result[1]


;delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/var_def,depth_def,lat_def,lon_def,var_new,lat_new,lon_new,depth_new/])

index=dimsizes(var_interp)

if(ivar.eq.0)then
	vars_3D = new((/dimsizes(var3D),index(0),index(1),index(2),index(3)/),typeof(var_interp))
	vars_3D!0 = "vars"
	vars_3D&vars = var3D
end if
vars_3D(ivar,:,:,:,:) = var_interp

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

fapend = (/"_"+season(is)+"_climo.nc", \
           "_"+season(is)+"_climo.nc", \
	   "_"+season(is)+"_climo.nc" /)

;;; for normal variables

;result_obs	= get_var_ocn_2D(fdir(0),fapend(0),casename(0),var2D(ivar),"obs")
;var_obs		= result_obs[0]
;lat_obs		= result_obs[1]
;lon_obs		= result_obs[2]

result_def	= get_var_ocn_2D(fdir(1),fapend(1),casename(1),var2D(ivar),"model")
var_def		= result_def[0]
lat_def		= result_def[1]
lon_def		= result_def[2]

result_new	= get_var_ocn_2D(fdir(2),fapend(2),casename(2),var2D(ivar),"model")
var_new		= result_new[0]
lat_new		= result_new[1]
lon_new		= result_new[2]


;result = regrid_diff(var_obs,lat_obs,lon_obs,\
;	             var_def,lat_def,lon_def,\
;		     var_new,lat_new,lon_new)
;var_interp = result[0]
;var_interp_diff = result[1]

result = regrid_diff_2D_2var(var_def,lat_def,lon_def,\
		             var_new,lat_new,lon_new)
var_interp = result[0]
var_interp_diff = result[1]


;delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/var_def,lat_def,lon_def,var_new,lat_new,lon_new/])

index=dimsizes(var_interp)

if(ivar.eq.0)then
	vars_2D = new((/dimsizes(var2D),index(0),index(1),index(2)/),typeof(var_interp))
	vars_2D!0 = "vars"
	vars_2D&vars = var2D
end if
vars_2D(ivar,:,:,:) = var_interp

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

vars_3D_mldmean		=	new((/dimsizes(var3D),ncase,nlat,nlon/),float)

; entrainment velocity, the vertical velocity at the mixed layer base
w_mld			=	new((/ncase,nlat,nlon/),float)

; the temperature difference between the mixed layer and just below the mixed layer
deltaT_mld		= 	new((/ncase,nlat,nlon/),float)

; for T, U, V, need to get the averaged values from sea surface to mixed layer depth.
; for W, should get the entrainment velocity, which is the vertical velocity at the mixed layer base.
do icase= 0,ncase-1
do ilat = 0,nlat-1
do ilon = 0,nlon-1

index_mld(icase,ilat,ilon) 		= 	closest_val(HMXL(icase,ilat,ilon),depth)	

if(.not.ismissing(index_mld(icase,ilat,ilon)))then

vars_3D_mldmean(:,icase,ilat,ilon) 	= dim_avg_n_Wrap(vars_3D(:,icase,:index_mld(icase,ilat,ilon),ilat,ilon),1)
w_mld(icase,ilat,ilon)			= vars_3D(3,icase,index_mld(icase,ilat,ilon),ilat,ilon)
deltaT_mld(icase,ilat,ilon)		= vars_3D(0,icase,index_mld(icase,ilat,ilon),ilat,ilon)-vars_3D(0,icase,index_mld(icase,ilat,ilon)+1,ilat,ilon)

else

vars_3D_mldmean(:,icase,ilat,ilon) 	=	vars_3D@_FillValue
w_mld(icase,ilat,ilon)			= 	vars_3D@_FillValue
deltaT_mld(icase,ilat,ilon)		=	vars_3D@_FillValue

end if

end do ; do ilon=
end do ; do ilat=
end do ; do icase=

printVarSummary(index_mld)
printVarSummary(vars_3D_mldmean)
printMinMax(vars_3D_mldmean(2,:,:,:),False)

printVarSummary(w_mld)
printMinMax(w_mld,False)

printVarSummary(deltaT_mld)
printMinMax(deltaT_mld,False)

;========================================================
; calculate entrainment term
;========================================================

entrain			= HMXL
entrain 		= w_mld*deltaT_mld/HMXL

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
;Temp_mldmean 		= 	vars_3D_mldmean(0,:,:,:)
Temp_mldmean 		= 	vars_3D(0,:,0,:,:)
gradLatLon  		= 	grad_latlon_cfd (Temp_mldmean, vars_3D&lat, vars_3D&lon, True, False)
                	               
;T_grad_lat  		= 	gradLatLon[0]    ; for clarity; explicitly extract variables from returned 'list' 
;T_grad_lon  		= 	gradLatLon[1]
T_grad_dy  		= 	gradLatLon[0] 
T_grad_dx  		= 	gradLatLon[1]


delete(gradLatLon)

; replace the generic attributes with more specific attributes

;T_grad_lon@long_name 	= 	"longitudinal gradient (derivative)"
;T_grad_lat@long_name 	= 	"latitudinal gradient (derivative)"
;T_grad_lat@units     	= 	"degC/m"
;T_grad_lon@units     	= 	"degC/m"

T_grad_dx@long_name 	= 	"zonal gradient"
T_grad_dy@long_name 	= 	"meridional gradient"
T_grad_dx@units     	= 	"degC/m"
T_grad_dy@units     	= 	"degC/m"

;###;////////////////////////////////////////////////////////////////////////////////
;###;////////////////////////////////////////////////////////////////////////////////
;###; QINYI: test gradient function to make sure that
;###; the extreme values of dTdx not coming from the grad_latlon_cfd function.
;###
;###TG			= 	Temp_mldmean
;###dimTG 			= 	dimsizes(TG)
;###
;###TGX_gradsg = new( dimTG, typeof(TG), getFillValue(TG) )  ; lon=>X
;###TGY_gradsg = new( dimTG, typeof(TG), getFillValue(TG) )  ; lat=>Y
;###gradsg(TG, TGX_gradsg, TGY_gradsg)     ; procedure for gaussian grids
;###
;###copy_VarCoords(TG, TGX_gradsg)         ; add meta data
;###copy_VarCoords(TG, TGY_gradsg) 
;###TGX_gradsg@long_name = "TGX: gradsg"
;###TGX_gradsg@units     = "K/s"
;###TGY_gradsg@long_name = "TGY: gradsg"
;###TGY_gradsg@units     = "K/s"
;###
;###print("")
;###printMinMax(TGY_gradsg,True )          ; unscaled
;###printMinMax(TGX_gradsg,False)
;###print("")
;###
;###;;;; NOTION HERE: "gradsg" does not handle missing values, and the input array must be on a global grid. If any missing values are encountered in a particular 2D input grid, then all of the values in the corresponding output grids will be set to the missing value defined by the output grids' _FillValue attributes.
;###;;; SO, I only can use the grad_latlon_cfd function here to get the SST gradient.
;###
;###;************************************************
;###; CENTERED FINITE DIFFERENCES:
;###; Use local function [grad_latlon_cfd]
;###; to compute zonal (X) and meridional (Y) gradients.
;###;************************************************
;###lat 			= 	TG&lat
;###lon			=	TG&lon
;###
;###gradLatLon = grad_latlon_cfd(TG, lat, lon, True, False)  
;###TGY_cfd    = gradLatLon[0]
;###TGX_cfd    = gradLatLon[1]
;###
;###TGY_cfd@long_name = "TGY: cfd"
;###TGY_cfd@units     = "K/s"
;###
;###TGX_cfd@long_name = "TGX: cfd"
;###TGX_cfd@units     = "K/s"
;###
;###printMinMax(TGY_cfd, True)             ; unscaled
;###printMinMax(TGX_cfd,False)
;###print("")
;###
;###;;; output some data
;###system("/bin/rm -f gradient-term.nc")
;###ncdf = addfile("gradient-term.nc" ,"c")
;###
;###ncdf->TGX_gradsg 	= TGX_gradsg
;###ncdf->TGY_gradsg	= TGY_gradsg
;###
;###ncdf->TGX_cfd		= TGX_cfd
;###ncdf->TGY_cfd		= TGY_cfd
;###
;###exit
;###;////////////////////////////////////////////////////////////////////////////////
;###;////////////////////////////////////////////////////////////////////////////////

; zonal advection: u_dTdx
;u 			= 	vars_3D(1,:,0,:,:)
u_mldmean		= 	vars_3D_mldmean(1,:,:,:)

;u_dTdx 			= 	T_grad_lat
;u_dTdx 			= 	u_mldmean*T_grad_lat
u_dTdx 			= 	T_grad_dx
u_dTdx 			= 	u_mldmean*T_grad_dx


;change unit from K/s to K/month
u_dTdx 			= 	u_dTdx*30*24*3600

u_dTdx@long_name	= 	"advection term: u_dTdx"
u_dTdx@units 		= 	"degC/month"

printVarSummary(u_dTdx)
printMinMax(u_dTdx,False)

; meridional advection: v_dTdy
;v 			= 	vars_3D(2,:,0,:,:)
v_mldmean		= 	vars_3D_mldmean(2,:,:,:)

;v_dTdy 			= 	T_grad_lon
;v_dTdy 			= 	v_mldmean*T_grad_lon
v_dTdy 			= 	T_grad_dy
v_dTdy 			= 	v_mldmean*T_grad_dy


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

sumall 			= netFlux-u_dTdx-v_dTdy-entrain
copy_VarCoords(netFlux,sumall)
sumall@longname 	= "sum of four terms"
sumall@units 		= "degC/month"

printVarSummary(sumall)
printMinMax(sumall(:,{-20:0},{180:240}),False)

;========================================================
; put all terms into one array to faciliate the plots
;========================================================
terms 			= (/"sumall","netFlux","u_dTdx","v_dTdy","entrain"/)

finalvar		= new((/5,ncase,nlat,nlon/),float)
finalvar!0		= "terms"
finalvar&terms		= terms

finalvar(0,:,:,:)	= netFlux
finalvar		= (/sumall,netFlux,(-1.)*u_dTdx,(-1.)*v_dTdy,(-1.)*entrain/)


finalvardiff 		= new((/5,nlat,nlon/),float)
finalvardiff!0		= "terms"
finalvardiff&terms	= terms

finalvardiff(0,:,:)	= netFlux(0,:,:)
finalvardiff		= (/sumall(1,:,:)-sumall(0,:,:),\
			    netFlux(1,:,:)-netFlux(0,:,:),(-1.)*u_dTdx(1,:,:)-(-1.)*u_dTdx(0,:,:),\
			    (-1.)*v_dTdy(1,:,:)-(-1.)*v_dTdy(0,:,:),(-1.)*entrain(1,:,:)-(-1.)*entrain(0,:,:)/)


; change units from degC/month to W/m2
factor			= 30.*24.*3600.
sumall_Wm2		= sumall*cp_sw*rho_sw*HMXL/factor
netFlux_Wm2		= netFlux*cp_sw*rho_sw*HMXL/factor
u_dTdx_Wm2		= (-1.)*u_dTdx*cp_sw*rho_sw*HMXL/factor
v_dTdy_Wm2		= (-1.)*v_dTdy*cp_sw*rho_sw*HMXL/factor
entrain_Wm2		= (-1.)*entrain*cp_sw*rho_sw*HMXL/factor

finalvar_Wm2		= finalvar
finalvar_Wm2		= (/sumall_Wm2,netFlux_Wm2,\
			    u_dTdx_Wm2,v_dTdy_Wm2,\
			    entrain_Wm2/)

finalvardiff_Wm2 	= finalvardiff
finalvardiff_Wm2	= (/sumall_Wm2(1,:,:)-sumall_Wm2(0,:,:),\
			    netFlux_Wm2(1,:,:)-netFlux_Wm2(0,:,:),u_dTdx_Wm2(1,:,:)-u_dTdx_Wm2(0,:,:),\
			    v_dTdy_Wm2(1,:,:)-v_dTdy_Wm2(0,:,:),entrain_Wm2(1,:,:)-entrain_Wm2(0,:,:)/)

printVarSummary(finalvar_Wm2)
printMinMax(finalvar_Wm2,False)

printVarSummary(finalvardiff_Wm2)
printMinMax(finalvardiff_Wm2,False)

;========================================================
; write into one NC file to check them
;========================================================
system("/bin/rm -f terms.nc")
ncdf = addfile("terms.nc" ,"c")

;filedimdef(ncdf,"time",-1,True)

ncdf->sumall 		= sumall
ncdf->netFlux 		= netFlux
ncdf->u_dTdx		= -1.*u_dTdx
ncdf->v_dTdy		= -1.*v_dTdy
ncdf->entrain		= -1.*entrain

ncdf->sumall_Wm2 	= sumall_Wm2
ncdf->netFlux_Wm2	= netFlux_Wm2
ncdf->u_dTdx_Wm2 	= u_dTdx_Wm2
ncdf->v_dTdy_Wm2 	= v_dTdy_Wm2
ncdf->entrain_Wm2	= entrain_Wm2


ncdf->T_grad_dx		= T_grad_dx
ncdf->T_grad_dy		= T_grad_dy
ncdf->Temp_mldmean	= Temp_mldmean
ncdf->u_mldmean		= u_mldmean
ncdf->v_mldmean		= v_mldmean
;ncdf->w_mld		= w_mld
ncdf->deltaT_mld	= deltaT_mld
ncdf->HMXL		= HMXL

exit

;*************************
; plotting parameters
;*************************
wks = gsn_open_wks("x11","./figure/ocn-heat-budget")

;gsn_define_colormap(wks,"WhBlGrYeRe")

;plot = new((/dimsizes(season),6/),graphic)

plot = new(5,graphic)

lat_min = -45.
lat_max = 50.
lon_min = 60.
lon_max = 290.
;lon_min = 0.
;lon_max = 360.


res = True
res@gsnFrame = False
res@gsnDraw = False
res@gsnSpreadColors = True

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
;res@gsnScalarContour = True

res@mpLandFillColor = "gray"
res@mpMinLonF = lon_min
res@mpMaxLonF = lon_max
res@mpMinLatF = lat_min
res@mpMaxLatF = lat_max

res@mpCenterLonF = 180

res@lbOrientation = "Vertical"
res@pmLabelBarOrthogonalPosF = -0.01

res@lbLabelStride    = 2                        ; every other label
res@tmXBTickSpacingF = 60.                      ; every 60 degrees 

res@mpGridAndLimbOn        = True               ; turn on grid lines
res@mpGridLineDashPattern = 5
res@mpGridLatSpacingF = 10.
res@mpGridLonSpacingF = 20.
res@mpGridMaskMode = "MaskNotOcean"

res@tmXBLabelFontHeightF  = 0.014
res@tmYLLabelFontHeightF  = 0.014

string_cli = (/"DEF","NEW"/)

resP = True
resP@gsnPanelRowSpec = True

;---------------------------------------------------------
; plot figures with units of Wm2
;---------------------------------------------------------

res@cnLevelSelectionMode = "ManualLevels" ; manual levels
res@cnMinLevelValF       = -100.
res@cnMaxLevelValF       = 100.
res@cnLevelSpacingF      = 10.

res@gsnRightString = "W/m2"

do icli = 0,1
res@gsnLeftString = string_cli(icli)
do ifin = 0,4
res@gsnCenterString = terms(ifin)
	plot(ifin) = gsn_csm_contour_map_ce(wks,finalvar_Wm2(ifin,icli,:,:),res) 
end do ; do ifin=
gsn_panel(wks,plot,(/1,2,2/),resP)
end do ; do icli=

; diff
res@cnLevelSelectionMode = "ManualLevels" ; manual levels
res@cnMinLevelValF       = -50.
res@cnMaxLevelValF       = 50.
res@cnLevelSpacingF      = 5.

string_diff = (/"NEW-DEF"/)
res@gsnLeftString = string_diff

do ifin = 0,4
res@gsnCenterString = terms(ifin)
	plot(ifin) = gsn_csm_contour_map_ce(wks,finalvardiff_Wm2(ifin,:,:),res) 
end do ; do ifin=
gsn_panel(wks,plot,(/1,2,2/),resP)


;---------------------------------------------------------
; plot figures with units of degC/month
;---------------------------------------------------------

res@cnLevelSelectionMode = "ManualLevels" ; manual levels
res@cnMinLevelValF       = -1.
res@cnMaxLevelValF       = 1.
res@cnLevelSpacingF      = 0.1

res@gsnRightString = "degC/month"

do icli = 0,1
res@gsnLeftString = string_cli(icli)
do ifin = 0,4
res@gsnCenterString = terms(ifin)
	plot(ifin) = gsn_csm_contour_map_ce(wks,finalvar(ifin,icli,:,:),res) 
end do ; do ifin=
gsn_panel(wks,plot,(/1,2,2/),resP)
end do ; do icli=

; diff
res@cnLevelSelectionMode = "ManualLevels" ; manual levels
res@cnMinLevelValF       = -0.5
res@cnMaxLevelValF       = 0.5
res@cnLevelSpacingF      = 0.1

string_diff = (/"NEW-DEF"/)
res@gsnLeftString = string_diff

do ifin = 0,4
res@gsnCenterString = terms(ifin)
	plot(ifin) = gsn_csm_contour_map_ce(wks,finalvardiff(ifin,:,:),res) 
end do ; do ifin=
gsn_panel(wks,plot,(/1,2,2/),resP)

exit


;------------------------------------------------------------------------------------------------------------

  ; add one box over eastern pacific that we care more about.
  dum1 = new(4,graphic)

  resp = True
  resp@gsLineColor = "blue"
  resp@gsLineThicknessF = 3.0

  ypts = (/30.0, 30.0, -30.0, -30.0, 30.0/)
  xpts = (/180.0, 220.0, 220.0, 180.0, 180.0/)

  ; draw each ine separately. Each line must contain two points.
  do i = 0,3
  	dum1(i) = gsn_add_polyline(wks,plot(is,5),xpts(i:i+1),ypts(i:i+1),resp)
  end do

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


gsn_panel(wks,plot(is,3:5),(/3,1/),False)

gsn_panel(wks,plot(is,:),(/2,3/),False)

delete(res)


;gsn_panel(wks,ndtooned(plot),(/5,6/),False)

delete(plot)
delete(fapend)

if(compare_type(icomp) .eq."amip-obs" .or. compare_type(icomp) .eq."cmip-obs")then
	break
end if

end do; do icase=

end do; do icomp=
                               
end
