load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************
;TEMP_levs		= (/-2,1,4,7,10,13,16,19,22,25,28,31,33/)
TEMP_levs		= (/7,10,13,16,19,22,25,28.0/)
;TEMP_dlevs		= (/-15,-12,-10,-8,-6,-2,0,2,6,8,10,12,15/)
TEMP_dlevs		= (/-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3/)
TEMP_dlevs_sp		= (/-3,-2.5,-2,-1.5,-1,1,1.5,2,2.5,3/)

;UVEL_levs		= (/-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12/)
UVEL_levs		= (/-12,-10,-8,-6,-4,-2,0,2/)
UVEL_dlevs		= (/-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,0,0.5,1.0,2.0,3.0,4.0,5.0/)
UVEL_dlevs_sp		= (/-5.0,-4.0,-3.0,-2.0,-1.0,1.0,2.0,3.0,4.0,5.0/)

;UVEL_dlevs		= (/-9.0,-7.0,-5.0,-3.0,-1.0,-0.5,0,0.5,1.0,3.0,5.0,7.0,9.0/)
;UVEL_dlevs_sp		= (/-9.0,-7.0,-5.0,-3.0,-1.0,1.0,3.0,5.0,7.0,9.0/)

;VVEL_levs		= (/-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0,1.0,2.0,3.0,4.0,5.0,6.0/)
VVEL_levs		= (/-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.,1.0/)
VVEL_dlevs		= (/-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5,1.8/)
VVEL_dlevs_sp		= (/-1.8,-1.5,-1.2,-0.9,-0.6,0.6,0.9,1.2,1.5,1.8/)

WVEL_levs		= (/-2,-1.5,-1.0,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,1.0,1.5,2.0/)
WVEL_dlevs		= (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6/)
WVEL_dlevs_sp		= (/-0.6,-0.5,-0.4,-0.3,-0.2,0.2,0.3,0.4,0.5,0.6/)

;WVEL_levs		= (/-0.2,-0.15,-1.0,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,1.0,1.5,2.0/)
;WVEL_dlevs		= (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6/)



;season=(/"ANN","DJF","MAM","JJA","SON"/)
season			=(/"ANN"/)
;var3D			=(/"TEMP","UVEL","VVEL","WVEL"/)
;var3D_levs		= (/TEMP_levs,UVEL_levs,VVEL_levs,WVEL_levs/)
;var3D_dlevs		= (/TEMP_dlevs,UVEL_dlevs,VVEL_dlevs,WVEL_dlevs/)
;var3D_units		= (/"degC","cm/s","cm/s","x10~S~-4~N~ cm s~S~-1~N~"/)

var3D			=(/"UVEL","VVEL","TEMP"/)
var3D_levs		= (/UVEL_levs,VVEL_levs,TEMP_levs/)
var3D_dlevs		= (/UVEL_dlevs,VVEL_dlevs,TEMP_dlevs/)
var3D_units		= (/"cm/s","cm/s","degC"/)
var3D_dlevs_sp		= (/UVEL_dlevs_sp,VVEL_dlevs_sp,TEMP_dlevs_sp/)


var2D			=(/"HMXL","SHF"/)
compare_type 		= (/"cmip-obs"/)
case			= (/"obs","def","new"/)

Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do is = 0,dimsizes(season)-1

ncdf = addfile("ocn-cross-section-w-obs-vars_3D_"+season(is)+".nc" ,"r")

dSizes                          = getfiledimsizes(ncdf)  

vars_3D         = new((/dimsizes(var3D),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)

do ivar=0,dimsizes(var3D)-1
        vars_3D(ivar,:,:,:,:)           =ncdf->$var3D(ivar)$ 

; change units from m/s to cm/s
	if(var3D(ivar).eq."UVEL".or.var3D(ivar).eq."VVEL")then
		vars_3D(ivar,:,:,:,:) 	= vars_3D(ivar,:,:,:,:)*1e2
	end if
end do

vars_3D!0		= "vars"

;========================================================
; Get the difference between two cases
;========================================================

vars_3D_diff			= vars_3D
vars_3D_diff(:,0,:,:,:)		= vars_3D(:,1,:,:,:)-vars_3D(:,0,:,:,:)
vars_3D_diff(:,1,:,:,:)		= vars_3D(:,2,:,:,:)-vars_3D(:,0,:,:,:)
vars_3D_diff(:,2,:,:,:)		= vars_3D(:,2,:,:,:)-vars_3D(:,1,:,:,:)

printVarSummary(vars_3D_diff)

;========================================================
; Prepare plot settings
;========================================================

latS_zm			= -20.
latE_zm			= -5.

lonS_zm			= 179.5
lonE_zm			= 260.5

VertSec_zm		= True


;========================================================
; zonal mean over specific zonal band
;========================================================

vars_3D_zm		= dim_avg_n_Wrap(vars_3D(:,:,:,{latS_zm:latE_zm},:),3)
printVarSummary(vars_3D_zm)

vars_3D_diff_zm		= dim_avg_n_Wrap(vars_3D_diff(:,:,:,{latS_zm:latE_zm},:),3)
printVarSummary(vars_3D_diff_zm)
printMinMax(vars_3D_diff_zm,True)


;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------

if(VertSec_zm)then

wks = gsn_open_wks("x11","./figure/ocn-cross-section-w-obs-zm_"+season(is)+"_revised-LonRange")

;gsn_define_colormap(wks,"BlueDarkRed18")
gsn_define_colormap(wks,"BlueWhiteOrangeRed")

plot 				= new((/dimsizes(case),dimsizes(var3D)/),graphic)
plotdiff 			= new((/dimsizes(case),dimsizes(var3D)/),graphic)

plot_ov				= new((/dimsizes(case),dimsizes(var3D)/),graphic)
plotdiff_ov 			= new((/dimsizes(case),dimsizes(var3D)/),graphic)

plot_ovv			= new((/2,dimsizes(var3D)/),graphic)


res 				= True

res@gsnFrame			= False
res@gsnDraw			= False

;res@vpWidthF			= 0.75
;res@vpHeightF			= 0.4

res@gsnSpreadColors 		= True
;res@gsnSpreadColorStart		= 20
;res@gsnSpreadColorEnd		= 237
res@gsnSpreadColorStart		= 40
res@gsnSpreadColorEnd		= 217


res@cnFillOn 			= True
res@cnLinesOn 			= False
res@cnLineLabelsOn 		= False

res@trYReverse			= True

res@lbOrientation 		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= -0.01

res@lbLabelStride    		= 2                        ; every other label
res@lbLabelFontHeightF		= 0.0325
res@tmXBTickSpacingF 		= 20

res@tmXBLabelFontHeightF  	= 0.0325
res@tmYLLabelFontHeightF  	= 0.0325

res@tiYAxisString		= "Depth (m)"
res@tiYAxisFontHeightF		= 0.0325

res@gsnStringFontHeightF	= 0.0325

res@tmYLMode			= "Explicit"
res@tmYLValues			= (/50,100,150,200,300,500,800/)
res@tmYLLabels			= (/50,100,150,200,300,500,800/)

res@gsnLeftString		= ""
res@gsnRightString		= ""
res@gsnCenterString		= ""

;-----------------------------------
; resource for later overlay
;-----------------------------------
resOV					= res
resOV@cnLinesOn				= True
resOV@cnFillOn				= False
resOV@gsnContourZeroLineThicknessF 	= 0
resOV@gsnContourNegLineDashPattern 	= 16
resOV@cnLineColor			= "black";"purple"
resOV@cnLineThicknessF			= 2.5

;---------------------------------------------------------
; plot figures of separate schemes
;---------------------------------------------------------

if(latE_zm.lt.0)then
res@gsnRightString 		= abs(latS_zm)+"S-"+abs(latE_zm)+"S"
else
res@gsnRightString 		= abs(latS_zm)+"S-"+abs(latE_zm)+"N"
end if


resOV2					= True	
resOV2@gsnFrame				= False
resOV2@gsnDraw				= False

resOV2@cnFillOn				= False
resOV2@cnLinesOn			= False
resOV2@cnLineLabelsOn			= False
resOV2@cnInfoLabelOn			= False

resOV2@gsnLeftString			= ""
resOV2@gsnRightString			= ""
resOV2@gsnCenterString			= ""

resOV2@cnSmoothingOn			= True
resOV2@cnSmoothingTensionF		= -10.0

;---------------------------------------------------------
; plot figures
;---------------------------------------------------------
;;; Pay attention: the first plotting here is the Shaded plot!!!

string_cli 			= (/"OBS","DEF","NEW"/)
string_diff 			= (/"DEF-OBS","NEW-OBS","NEW-DEF"/)

do ivar = 0,dimsizes(var3D)-1

	res@gsnCenterString 		= var3D(ivar)
	res@gsnRightString		= var3D_units(ivar)


;;; Pay attention: the overlaying plot here is the contour plot!!!

	do icase = 0,2

		; climatology  -- shaded
		res@cnFillOn			= True
		res@cnLinesOn			= False
		res@cnLineLabelsOn		= False
		res@gsnLeftString		= string_cli(icase)

		res@cnLevelSelectionMode        = "ExplicitLevels"
		res@cnLevels                    = var3D_levs(ivar,:)
		plot(icase,ivar)		= gsn_csm_contour(wks,vars_3D_zm(ivar,icase,{:1000},{lonS_zm:lonE_zm}),res) 
		delete(res@cnLevels)

		; overlaying climatology -- contour
		resOV@cnFillOn			= False
		resOV@cnLinesOn			= True
		resOV@cnLineLabelsOn		= False
		resOV@cnLineLabelFontHeightF	= 0.03
		resOV@cnLineLabelPlacementMode	= "Constant"
		resOV@cnLineDashSegLenF		= 0.1
		resOV@cnInfoLabelOn		= True
		resOV@cnInfoLabelFontHeightF	= 0.030

		resOV@cnLevelSelectionMode	= "ExplicitLevels"
		resOV@cnLevels			= var3D_levs(ivar,:)
		p				= 0.50
		q				= 0.25
;		tmp				= smth9(vars_3D_zm(ivar,icase,{:1000},{lonS_zm:lonE_zm}),p,q,False)
		tmp				= vars_3D_zm(ivar,icase,{:1000},{lonS_zm:lonE_zm})
		plot_ov(icase,ivar)		= gsn_csm_contour(wks,tmp,resOV) 
;		plot_ov(icase,ivar)		= gsn_csm_contour(wks,vars_3D_zm(ivar,icase,{:1000},{lonS_zm:lonE_zm}),resOV) 

		; these two plots are just for overlaying on DEF-OBS, and NEW-OBS.
		; this is temporary. if finding other good methods, pls update it.
;		tmp				= smth9(vars_3D_zm(ivar,0,{:1000},{lonS_zm:lonE_zm}),p,q,False)
		tmp				= vars_3D_zm(ivar,0,{:1000},{lonS_zm:lonE_zm})
		plot_ovv(0,ivar)		= gsn_csm_contour(wks,tmp,resOV)
		plot_ovv(1,ivar)		= gsn_csm_contour(wks,tmp,resOV)
		delete(resOV@cnLevels)

		; differences -- shaded
		; overlaying differences -- contour

		res@gsnLeftString 		= string_diff(icase)
		res@cnLevelSelectionMode	= "ExplicitLevels"
		res@cnLevels			= var3D_dlevs(ivar,:)
		p				= 0.50
		q				= 0.25
;		tmp				= smth9(vars_3D_diff_zm(ivar,icase,{:1000},{lonS_zm:lonE_zm}),p,q,False)
		tmp				= vars_3D_diff_zm(ivar,icase,{:1000},{lonS_zm:lonE_zm})
		plotdiff(icase,ivar)		= gsn_csm_contour(wks,tmp,res) 
		delete(res@cnLevels)

		resOV2@cnLevelSelectionMode	= "ExplicitLevels"
		resOV2@cnLevels			= var3D_dlevs_sp(ivar,:)
		plotdiff_ov(icase,ivar)		= gsn_csm_contour(wks,tmp,resOV2) 
		overlay(plotdiff(icase,ivar),plotdiff_ov(icase,ivar))
		delete(resOV2@cnLevels)

		if(icase.ne.2)then
			overlay(plotdiff(icase,ivar),plot_ovv(icase,ivar))
			;---Remove plotA from map1 before overlaying it again
			; NhlRemoveOverlay(plotdiff(icase,ivar),plot_ov(0,ivar),False) 
		end if
	end do ; do icase=

end do ; do ivar=

resP1 				= True
resP1@gsnPanelXWhiteSpacePercent 	= 5.
resP1@gsnPanelYWhiteSpacePercent 	= 5.
resP1@amJust                     = "BottomLeft"
resP1@gsnPanelFigureStrings      = (/"a)","b)","c)","d)","e)","f)","g)","h)","i)"/)

; U,V,T
plotall_1		= (/plotdiff(:,0),plotdiff(:,1),plotdiff(:,2)/)
; T,U,V
plotall_2		= (/plotdiff(:,2),plotdiff(:,0),plotdiff(:,1)/)

gsn_panel(wks,ndtooned(plotall_1),(/3,3/),resP1)
gsn_panel(wks,ndtooned(plotall_2),(/3,3/),resP1)

delete(plotall_1)
delete(plotall_2)

delete(res)
delete(plot)
delete(wks)

end if ; VertSec_zm

end do ; do is=
                               
end
