load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************
TEMP_levs		= (/4,6,8,10,12,14,16,18,20,25,26,28,30/)
;TEMP_dlevs		= (/-15,-12,-10,-8,-6,-2,0,2,6,8,10,12,15/)
TEMP_dlevs		= (/-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3/)

UVEL_levs		= (/-11,-9,-7,-5,-3,-1,0,1,3,5,7,9,11/)
UVEL_dlevs		= (/-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,0,0.5,1.0,2.0,3.0,4.0,5.0/)

VVEL_levs		= (/-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,0,0.5,1.0,2.0,3.0,4.0,5.0/)
VVEL_dlevs		= (/-2,-1.5,-1.2,-1.0,-0.5,-0.3,0,0.3,0.5,1.0,1.2,1.5,2.0/)

WVEL_levs		= (/-2,-1.5,-1.0,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,1.0,1.5,2.0/)
WVEL_dlevs		= (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6/)

;WVEL_levs		= (/-0.2,-0.15,-1.0,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,1.0,1.5,2.0/)
;WVEL_dlevs		= (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6/)



;season=(/"ANN","DJF","MAM","JJA","SON"/)
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

vars_3D         = new((/dimsizes(var3D),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)

do ivar=0,dimsizes(var3D)-1
        vars_3D(ivar,:,:,:,:)           =ncdf->$var3D(ivar)$ 
end do

vars_3D!0		= "vars"

printVarSummary(vars_3D)

;change the unit of w from cm/s to 10^-3 cm/s
;x10~S~-3~N~ cm s~S~-1~N~
vars_3D(3,:,:,:,:)	= vars_3D(3,:,:,:,:)*1e4

;========================================================
; Get the difference between two cases
;========================================================

vars_3D_diff		= vars_3D(:,1,:,:,:)-vars_3D(:,0,:,:,:)

copy_VarCoords_1(vars_3D(vars|:,lev|:,lat|:,lon|:,case|:),vars_3D_diff)
printVarSummary(vars_3D_diff)


;========================================================
; Prepare plot settings
;========================================================

latS_zm			= -15.
latE_zm			= -5.

lonS_zm			= 160.
lonE_zm			= 280.

VertSec_zm		= True


;========================================================
; zonal mean over specific zonal band
;========================================================

vars_3D_zm		= dim_avg_n_Wrap(vars_3D(:,:,:,{latS_zm:latE_zm},:),3)
printVarSummary(vars_3D_zm)

vars_3D_diff_zm		= dim_avg_n_Wrap(vars_3D_diff(:,:,{latS_zm:latE_zm},:),2)
printVarSummary(vars_3D_diff_zm)
printMinMax(vars_3D_diff_zm,True)




;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------

if(VertSec_zm)then

wks = gsn_open_wks("x11","./figure/ocn-cross-section-zm_"+season)

plot 				= new((/dimsizes(case),dimsizes(var3D)/),graphic)
plotdiff 			= new((/dimsizes(var3D)/),graphic)

plot_ov				= new((/dimsizes(case),dimsizes(var3D)/),graphic)
plotdiff_ov 			= new((/dimsizes(var3D)/),graphic)


res 				= True

res@gsnFrame			= False
res@gsnDraw			= False

;res@vpWidthF			= 0.75
;res@vpHeightF			= 0.4

res@gsnSpreadColors 		= True

res@cnFillOn 			= True
res@cnLinesOn 			= True
res@cnLineLabelsOn 		= False

res@trYReverse			= True


res@lbOrientation 		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= -0.01

res@lbLabelStride    		= 2                        ; every other label
res@tmXBTickSpacingF 		= 30

res@tmXBLabelFontHeightF  	= 0.025
res@tmYLLabelFontHeightF  	= 0.025

res@tiYAxisString		= "Depth (m)"

res@gsnStringFontHeightF	= 0.025


string_cli 			= (/"DEF","NEW"/)

resP 				= True
resP@amJust                     = "BottomLeft"
resP@gsnPanelFigureStrings      = (/"a)","b)","c)"/)


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

if(latE_zm.lt.0)then
res@gsnRightString 		= abs(latS_zm)+"S-"+abs(latE_zm)+"S"
else
res@gsnRightString 		= abs(latS_zm)+"S-"+abs(latE_zm)+"N"
end if


do ivar = 0,dimsizes(var3D)-1

if(ivar.ne.0)then
	res@tiYAxisString			= ""
	resOV@tiYAxisString			= ""
end if

	res@gsnCenterString 		= var3D(ivar)

	res@cnLevelSelectionMode        = "ExplicitLevels"
	res@cnLevels                    = var3D_levs(ivar,:)

	res@gsnRightString		= var3D_units(ivar)

	do icase = 0,1
		res@cnLinesOn			= True
		res@gsnLeftString 		= string_cli(icase)
		res@cnFillOn			= True
		plot(icase,ivar)	= gsn_csm_contour(wks,vars_3D_zm(ivar,icase,{:1000},{lonS_zm:lonE_zm}),res) 

		; draw climatology as lines for overlaying later.
		plot_ov(icase,ivar)	= gsn_csm_contour(wks,vars_3D_zm(ivar,icase,{:1000},{lonS_zm:lonE_zm}),resOV) 
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

if(ivar.ne.0)then
	res@tiYAxisString			= ""
else
	res@tiYAxisString			= "Depth(m)"
end if

	res@cnLevelSelectionMode        = "ExplicitLevels"
	res@cnLevels                    = var3D_dlevs(ivar,:)

	res@gsnCenterString 		= var3D(ivar)
	res@cnLinesOn			= True
	res@cnFillOn			= True

	res@gsnRightString		= var3D_units(ivar)

	plotdiff(ivar)		= gsn_csm_contour(wks,vars_3D_diff_zm(ivar,{:1000},{lonS_zm:lonE_zm}),res) 

end do ; do ivar=

;gsn_panel(wks,plotdiff,(/2,2/),resP)

;gsn_panel(wks,plotdiff(:2),(/1,3/),resP)

plotall				= (/plot(0,:2),plotdiff(:2)/)

;gsn_panel(wks,ndtooned(plotall),(/2,3/),resP)

;---------------------------------------------------------
; Add overlay plots of mean climate and their difference
;---------------------------------------------------------

resOV2				= res

do ivar=0,dimsizes(var3D)-1

if(ivar.ne.0)then
	resOV2@tiYAxisString			= ""
else
	resOV2@tiYAxisString			= "Depth(m)"
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

	plotdiff_ov(ivar)			= gsn_csm_contour(wks,vars_3D_diff_zm(ivar,{:1000},{lonS_zm:lonE_zm}),resOV2) 

	overlay(plotdiff_ov(ivar),plot_ov(0,ivar))
end do
;gsn_panel(wks,plot_ov(0,:2),(/1,3/),resP)
gsn_panel(wks,plotdiff_ov(:2),(/1,3/),resP)


delete(res)
delete(plot)
delete(wks)

end if ; VertSec_zm

end do ; do is=
                               
end
