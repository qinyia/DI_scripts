load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************

TEMP_levs		= (/4.,6,8,10,12,14,16,18,20,25,26,28,30/)
TEMP_dlevs		= (/-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3/)

UVEL_levs		= (/-80.,-70,-60,-40,-30,-20,0,20,30,40,60,70,80/)
UVEL_dlevs		= (/-15.0,-12.0,-9.0,-6.0,-3.0,-1,0,1,3.0,6.0,9.0,12.0,15.0/)

VVEL_levs		= (/-20.,-15,-10,-6,-4,-2,0,2,4,6,10,15,20/)
VVEL_dlevs		= (/-2,-1.5,-1.2,-1.0,-0.5,-0.3,0,0.3,0.5,1.0,1.2,1.5,2.0/)

WVEL_levs		= (/-20.,-15,-10,-6,-4,-2,0,2,4,6,10,15,20/)
WVEL_dlevs		= (/-6.0,-5.0,-4,-3,-2,-1,0,1,2,3,4,5,6/)



season			=(/"ANN"/)

var3D			=(/"TEMP","UVEL","VVEL","WVEL"/)
var3D_levs		= (/TEMP_levs,UVEL_levs,VVEL_levs,WVEL_levs/)
var3D_dlevs		= (/TEMP_dlevs,UVEL_dlevs,VVEL_dlevs,WVEL_dlevs/)
var3D_units		= (/"degC","cm/s","cm/s","x10~S~-4~N~ cm s~S~-1~N~"/)


var2D			=(/"HMXL","SHF"/)
compare_type 		= (/"cmip-obs"/)
case			= (/"def","new"/)

Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do is = 0,dimsizes(season)-1

ncdf = addfile("ocn-cross-section-vars_3D_"+season(is)+".nc" ,"r")

dSizes                          = getfiledimsizes(ncdf)  
print(dSizes)
vars_3D         = new((/dimsizes(var3D),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)

do ivar=0,dimsizes(var3D)-1
        vars_3D(ivar,:,:,:,:)           =ncdf->$var3D(ivar)$ 
end do

vars_3D!0		= "vars"

; w: cm/s to 10^-4 cm/s
vars_3D(3,:,:,:,:)		= vars_3D(3,:,:,:,:)*1e4

printVarSummary(vars_3D)

;========================================================
; Get the difference between two cases
;========================================================

vars_3D_diff		= vars_3D(:,1,:,:,:)-vars_3D(:,0,:,:,:)

copy_VarCoords_1(vars_3D(vars|:,lev|:,lat|:,lon|:,case|:),vars_3D_diff)
printVarSummary(vars_3D_diff)


;========================================================
; Prepare plot settings
;========================================================

latS_mm			= -20.
latE_mm			=  20.

lonS_mm			= 180.
lonE_mm			= 210.

VertSec_mm		= True


;========================================================
; meridional mean over specific meridional band
;========================================================

vars_3D_mm		= dim_avg_n_Wrap(vars_3D(:,:,:,:,{lonS_mm:lonE_mm}),4)
printVarSummary(vars_3D_mm)

vars_3D_diff_mm		= dim_avg_n_Wrap(vars_3D_diff(:,:,:,{lonS_mm:lonE_mm}),3)
printVarSummary(vars_3D_diff_mm)

do ivar=0,dimsizes(var3D)-1
printMinMax(vars_3D_mm(ivar,:,:,:),True)
printMinMax(vars_3D_diff_mm(ivar,:,:),True)
end do

;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------
 
if(VertSec_mm)then

wks = gsn_open_wks("x11","./figure/ocn-cross-section-mm")

plot 				= new((/dimsizes(case),dimsizes(var3D)/),graphic)
plotdiff 			= new((/dimsizes(var3D)/),graphic)

plot_ov				= new((/dimsizes(case),dimsizes(var3D)/),graphic)
plotdiff_ov 			= new((/dimsizes(var3D)/),graphic)


res 				= True

res@gsnFrame			= False
res@gsnDraw			= False

;res@vpWidthF			= 0.75
;res@vpHeightF			= 0.40

res@gsnSpreadColors 		= True

res@cnFillOn 			= True
res@cnLinesOn 			= True
res@cnLineLabelsOn 		= False

res@trYReverse			= True


res@lbOrientation 		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= -0.01

res@lbLabelStride    		= 2                        ; every other label
res@tmXBTickSpacingF 		= 15.


res@tmXBLabelFontHeightF  	= 0.025
res@tmYLLabelFontHeightF  	= 0.025

res@gsnStringFontHeightF	= 0.025

string_cli 			= (/"DEF","NEW"/)

resP 				= True
resP@amJust                     = "BottomLeft"
resP@gsnPanelFigureStrings      = (/"d)","e)","f)"/)

;-----------------------------------
; resource for later overlay
;-----------------------------------
resOV					= res
resOV@cnLinesOn				= True
resOV@cnFillOn				= False
resOV@gsnContourZeroLineThicknessF 	= 0
resOV@gsnContourNegLineDashPattern 	= 16
resOV@cnLineColor			= "white"
resOV@cnLineThicknessF			= 3.0



;---------------------------------------------------------
; plot figures of separate schemes
;---------------------------------------------------------

res@gsnRightString 		= lonS_mm+"E-"+lonE_mm+"E"


do ivar = 0,dimsizes(var3D)-1

if(ivar.eq.0)then
	res@tiYAxisString		= "Depth(m)"
	resOV@tiYAxisString		= "Depth(m)"
else
	res@tiYAxisString		= ""
	resOV@tiYAxisString		= ""
end if

	res@gsnCenterString 		= var3D(ivar)

	res@cnLevelSelectionMode        = "ExplicitLevels"
	res@cnLevels                    = var3D_levs(ivar,:)

	do icase = 0,1
		res@gsnLeftString 		= string_cli(icase)
		plot(icase,ivar)	= gsn_csm_contour(wks,vars_3D_mm(ivar,icase,{:1000},{latS_mm:latE_mm}),res) 

		; draw climatology as lines for overlaying later.
		plot_ov(icase,ivar)	= gsn_csm_contour(wks,vars_3D_mm(ivar,icase,{:1000},{latS_mm:latE_mm}),resOV) 

	end do ; do icase=
;	gsn_panel(wks,plot(:,ivar),(/1,2/),resP)
end do ; do ivar=


delete(ivar)
delete(icase)


;---------------------------------------------------------
; plot difference figures
;---------------------------------------------------------

string_diff 			= (/"NEW-DEF"/)
res@gsnLeftString 		= string_diff

do ivar = 0,dimsizes(var3D)-1

if(ivar.eq.0)then
	res@tiYAxisString		= "Depth(m)"
else
	res@tiYAxisString		= ""
end if

	res@cnLevelSelectionMode	= "ExplicitLevels"
	res@cnLevels			= var3D_dlevs(ivar,:)

	res@gsnCenterString 		= var3D(ivar)

	plotdiff(ivar)		= gsn_csm_contour(wks,vars_3D_diff_mm(ivar,{:1000},{latS_mm:latE_mm}),res) 

end do ; do ivar=

;gsn_panel(wks,plotdiff,(/2,2/),resP)

;---------------------------------------------------------
; Add overlay plots of mean climate and their difference
;---------------------------------------------------------

resOV2				= res

do ivar=0,dimsizes(var3D)-1

if(ivar.eq.0)then
	resOV2@tiYAxisString		= "Depth(m)"
else
	resOV2@tiYAxisString		= ""
end if




	resOV2@cnLevelSelectionMode        	= "ExplicitLevels"
	resOV2@cnLevels                    	= var3D_dlevs(ivar,:)

	resOV2@cnFillOn				= True
	resOV2@cnLineLabelsOn			= False
	resOV2@cnLinesOn			= True

	resOV2@gsnCenterString			= ""
	resOV2@gsnRightString			= ""

	resOV2@gsnLeftString			= "NEW-DEF"
	resOV2@gsnCenterString			= var3D(ivar)
	resOV2@gsnRightString			= var3D_units(ivar)

	plotdiff_ov(ivar)			= gsn_csm_contour(wks,vars_3D_diff_mm(ivar,{:1000},{latS_mm:latE_mm}),resOV2) 

	overlay(plotdiff_ov(ivar),plot_ov(0,ivar))
end do
;gsn_panel(wks,plot_ov(0,:2),(/1,3/),resP)
gsn_panel(wks,plotdiff_ov(:2),(/1,3/),resP)


delete(res)
delete(plot)
delete(wks)



end if ; VertSec_mm

end do ; do is=
 
                               
end
