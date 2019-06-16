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

DataProcess		= False

;season=(/"ANN","DJF","MAM","JJA","SON"/)

season			=(/"ANN"/)
;var3D			=(/"TEMP","UVEL","VVEL","WVEL"/)
var3D			=(/"TEMP"/)

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

if(DataProcess)then

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

;========================================================
; calculate horizontal advection: zonal and meridional
;========================================================
; get latitudinal and longitudinal gradient
; here, leave one small problem: how to directly use "TEMP" as the index?
; temperature is the first variable in var3D, so its index is 0.
Temp_mldmean 		= 	vars_3D(0,:,0,:,:)

gradLatLon  		= 	grad_latlon_cfd (Temp_mldmean, vars_3D&lat, vars_3D&lon, True, False)
                	               
T_grad_dy  		= 	gradLatLon[0] 
T_grad_dx  		= 	gradLatLon[1]

delete(gradLatLon)

T_grad_dx@long_name 	= 	"zonal gradient"
T_grad_dy@long_name 	= 	"meridional gradient"
T_grad_dx@units     	= 	"degC/m"
T_grad_dy@units     	= 	"degC/m"

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

;========================================================
; write into one NC file to check them
;========================================================
system("/bin/rm -f gradient"+"_"+season+".nc")
ncdf = addfile("gradient"+"_"+season+".nc" ,"c")

ncdf->T_grad_dx		= T_grad_dx
ncdf->T_grad_dy		= T_grad_dy

ncdf->T_grad_dx_diff	= T_grad_dx_diff
ncdf->T_grad_dy_diff	= T_grad_dy_diff

end if

;========================================================
; Read data and Calculate zonal mean
;========================================================

ncin = addfile("gradient"+"_"+season+".nc" ,"w")

dTdx			= ncin->T_grad_dx

printVarSummary(dTdx)

dTdx			= dTdx*1e6	; K/m -> 10^-3 K/km

latS					= -20.
latE					= -5.
lonS					= 160.
lonE					= 260.

lat_here                                = dTdx&lat
rad                                     = 4.0*atan(1.0)/180.
clat                                    = cos(lat_here*rad)
clat!0                                  = "lat"
clat&lat                                = lat_here

dTdx_mm					= dim_avg_wgt_n_Wrap(dTdx(:,{latS:latE},:),clat({latS:latE}),1,1)

printVarSummary(dTdx_mm)

;========================================================
; Parameters Setting...
;========================================================
LonValuePlot		= True
PlotType		= "x11"


if(LonValuePlot)then
;---------------------------------------------------------------
; Draw Lon-Value Plots: separate schemes
;---------------------------------------------------------------
wks = gsn_open_wks(PlotType,"./figure/SST-gradient_"+season)

plot				=new(2,graphic)

resL 				= True

resL@vpWidthF			= 0.7
resL@vpHeightF			= 0.4

resL@gsnFrame		 	= False
resL@gsnDraw			= False

resL@xyLineThicknessF 		= 3.0

resL@tmXBTickSpacingF           = 15.

resL@gsnYRefLine       		= 0.0
resL@gsnYRefLineColor		= "grey" 
resL@gsnYRefLineThicknessF	= 3.0

resL@tiYAxisString		= "SST gradient(10~S~-3~N~ K km~S~-1~N~)"

;resL@xyDashPatterns		= (/0,0,0,0,0/)
resL@xyLineColors		= (/"red","blue"/)

tmp				= dTdx_mm(:,{lonS:lonE})

plot(0)			 	= gsn_csm_xy(wks,tmp&lon,tmp,resL)
;----------------------------------------
; add legends
;----------------------------------------

legends                         = (/"DEF","NEW"/)

lgres                           = True
lgres@lgLineColors            	= resL@xyLineColors
lgres@lgLineThicknessF		= resL@xyLineThicknessF
;lgres@lgDashIndexes		= resL@xyDashPatterns
lgres@lgLabelFontHeightF        = 0.08
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.1
lgres@lgPerimOn                 = False

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres                           = True
amres@amJust			= "BottomLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = 0.50

annoid1                         = gsn_add_annotation(plot(0),lbid1,amres)

gsn_panel(wks,plot(0),(/1,1/),False)

delete(resL)
delete(plot)
delete(wks)

end if ; LonValuePlot



end do; do icomp=
                               
end
