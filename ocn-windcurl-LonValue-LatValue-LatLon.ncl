load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************

;season=(/"ANN","DJF","MAM","JJA","SON"/)
season			=(/"ANN"/)
var3D			=(/"TEMP","UVEL","VVEL","WVEL"/)
var2D			=(/"TAUX","TAUY"/)
compare_type 		= (/"cmip-obs"/)
case			= (/"def","new"/)

Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp=0,dimsizes(compare_type)-1
do is=0,dimsizes(season)-1

ncdf = addfile("ocn-var3D_"+compare_type(icomp)+"_"+season(is)+"_wo_obs.nc","r")

dSizes                          = getfiledimsizes(ncdf)  

vars_LevLatLon         = new((/dimsizes(var3D),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)

do ivar=0,dimsizes(var3D)-1
        vars_LevLatLon(ivar,:,:,:,:)           =ncdf->$var3D(ivar)$ 
end do

vars_LevLatLon!0		= "vars"

delete(dSizes)

printVarSummary(vars_LevLatLon)
printMinMax(vars_LevLatLon,True)

; for 2D varialbes
ncdf = addfile("ocn-var2D_"+compare_type(icomp)+"_"+season(is)+"_wo_obs.nc","r")

dSizes                          = getfiledimsizes(ncdf)  

vars_LatLon         = new((/dimsizes(var2D),dSizes(0),dSizes(1),dSizes(2)/),float)

do ivar=0,dimsizes(var2D)-1
        vars_LatLon(ivar,:,:,:)           =ncdf->$var2D(ivar)$ 
end do

vars_LatLon!0		= "vars"

printVarSummary(vars_LatLon)
printMinMax(vars_LatLon,True)

;========================================================
; Prepare plot settings
;========================================================

latS			= -30.
latE			= 30

lonS			= 120.
lonE			= 280.

levS			= 50.

lev			= (/5/)


WindCurlPlot		= False
LonValuePlot		= False
LatValuePlot		= False

LatLonPlot		= True

PlotType		= "x11"


;========================================================
; Get the difference between two cases
;========================================================

vars_LevLatLon_diff		= vars_LevLatLon(:,1,:,:,:)-vars_LevLatLon(:,0,:,:,:)
copy_VarCoords_1(vars_LevLatLon(vars|:,lev|:,lat|:,lon|:,case|:),vars_LevLatLon_diff)
printVarSummary(vars_LevLatLon_diff)

vars_LatLon_diff		= vars_LatLon(:,1,:,:)-vars_LatLon(:,0,:,:)
copy_VarCoords_1(vars_LatLon(vars|:,lat|:,lon|:,case|:),vars_LatLon_diff)
printVarSummary(vars_LatLon_diff)

;========================================================
; Get the wind stress curl
;========================================================
taux				= vars_LatLon(0,:,:,:)
tauy				= vars_LatLon(1,:,:,:)
uv2vr				= uv2vr_cfd(taux,tauy,taux&lat,taux&lon,2)
uv2vr				= uv2vr*1e7	; N/m2 -> 10^-7 N/m2
copy_VarCoords(taux,uv2vr)

uv2vr_diff			= uv2vr(1,:,:)-uv2vr(0,:,:)
uv2vr_diff			= uv2vr_diff
copy_VarCoords_1(uv2vr(lat|:,lon|:,case|:),uv2vr_diff)
printMinMax(uv2vr_diff,True)



;========================================================
; Get the meridional mean over specific latitude ranges--> Lon-value plots
;========================================================

vars_LatLon_mm			= dim_avg_n_Wrap(vars_LatLon(:,:,{latS:latE},:),2)
printVarSummary(vars_LatLon_mm)

vars_LatLon_diff_mm		= dim_avg_n_Wrap(vars_LatLon_diff(:,{latS:latE},:),1)
printVarSummary(vars_LatLon_diff_mm)
printMinMax(vars_LatLon_diff_mm,True)

;========================================================
; Get the zonal mean over specific longitude ranges--> Lat-Value plots
;========================================================
vars_LatLon_zm			= dim_avg_n_Wrap(vars_LatLon(:,:,:,{lonS:lonE}),3)
printVarSummary(vars_LatLon_zm)

vars_LatLon_diff_zm		= dim_avg_n_Wrap(vars_LatLon_diff(:,:,{lonS:lonE}),2)
printVarSummary(vars_LatLon_diff_zm)
printMinMax(vars_LatLon_diff_zm,True)


;--------------------------------------------------------------------------------------------------------
; Plotting starts here...
;--------------------------------------------------------------------------------------------------------

if(WindCurlPlot)then

;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------
wks = gsn_open_wks(PlotType,"./figure/ocn-lat-lon-windcurl-overlay-windspeed")

plot 				= new(3,graphic)


res 				= True

res@gsnFrame			= False
res@gsnDraw			= False

res@gsnSpreadColors 		= True

res@cnFillOn 			= True
res@cnLinesOn 			= True
res@cnLineLabelsOn 		= False

res@lbOrientation 		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= -0.01

res@lbLabelStride    		= 2                        ; every other label
;res@tmXBTickSpacingF 		= 60.                      ; every 60 degrees 
;res@tmXBTickSpacingF 		= 15.


res@tmXBLabelFontHeightF  	= 0.014
res@tmYLLabelFontHeightF  	= 0.014

res@mpCenterLonF		= 180.
res@mpMinLatF			= latS
res@mpMaxLatF			= latE
res@mpMinLonF			= lonS
res@mpMaxLonF			= lonE

res@gsnAddCyclic		= False

res@gsnScalarContour     	= True               ; contours desired

res@vcRefMagnitudeF		= 1.0             ; define vector ref mag
res@vcRefLengthF		= 0.045           ; define length of vec ref
res@vcRefAnnoOrthogonalPosF	= -1.0            ; move ref vector
res@vcRefAnnoArrowLineColor	= "black"         ; change ref vector color
res@vcRefAnnoArrowUseVecColor	= False           ; don't use vec color for ref

;res@vcGlyphStyle		= "CurlyVector"     ; turn on curly vectors
res@vcLineArrowColor		= "white"           ; change vector color
res@vcLineArrowThicknessF	= 4.0               ; change vector thickness
res@vcVectorDrawOrder		= "PostDraw"        ; draw vectors last
res@vcMinDistanceF		= 0.025

res@vcMinMagnitudeF		= 0.5		; NOTE: this cannot work when "CurlyVector" is turned on.

res@gsnRightString		= "Wind stress curl(shaded) & wind speed (vector)"

string_cli 			= (/"DEF","NEW"/)

resP 				= True
;resP@gsnPanelRowSpec 		= True


;---------------------------------------------------------
; plot figures of separate schemes
;---------------------------------------------------------
res@cnLevelSelectionMode	= "ManualLevels"
res@cnMinLevelValF		= -20.
res@cnMaxLevelValF		=  20.
res@cnLevelSpacingF		=   4.

do icase = 0,1
	res@gsnLeftString 		= string_cli(icase)
	taux_tmp			= vars_LatLon(0,icase,{latS:latE},{lonS:lonE})
	tauy_tmp			= vars_LatLon(1,icase,{latS:latE},{lonS:lonE})
	uv2vr_tmp			= uv2vr(icase,{latS:latE},{lonS:lonE})

	plot(icase)	= gsn_csm_vector_scalar_map_ce(wks,taux_tmp,tauy_tmp,uv2vr_tmp,res)

end do ; do icase=

;gsn_panel(wks,plot(0:1),(/2,1/),resP)

delete(taux_tmp)
delete(tauy_tmp)
delete(uv2vr_tmp)

;---------------------------------------------------------
; plot difference figures
;---------------------------------------------------------
res@cnLevelSelectionMode	= "ManualLevels"
res@cnMinLevelValF		= -8.
res@cnMaxLevelValF		= 8.
res@cnLevelSpacingF		= 1.


string_diff 			= (/"NEW-DEF"/)
res@gsnLeftString 		= string_diff

res@vcMinMagnitudeF		= 0.20		; NOTE: this cannot work when "CurlyVector" is turned on.
res@vcRefMagnitudeF		= 0.5             ; define vector ref mag

taux_tmp			= vars_LatLon_diff(0,{latS:latE},{lonS:lonE})
tauy_tmp			= vars_LatLon_diff(1,{latS:latE},{lonS:lonE})
uv2vr_tmp			= uv2vr_diff({latS:latE},{lonS:lonE})

plot(2)	= gsn_csm_vector_scalar_map_ce(wks,taux_tmp,tauy_tmp,uv2vr_tmp,res)

gsn_panel(wks,plot,(/3,1/),resP)

delete(res)
delete(plot)
delete(wks)


end if ; WindCurlPlot



;---------------------------------------------------------------
; Draw Lon-Value Plots
;---------------------------------------------------------------
if(LonValuePlot)then
;---------------------------------------------------------------
; Draw Lon-Value Plots: separate schemes
;---------------------------------------------------------------
wks = gsn_open_wks(PlotType,"./figure/ocn-LonValue_"+latS+"_"+latE)

plot				=new(2,graphic)

res 				= True

res@gsnFrame		 	= False
res@gsnDraw			= False

res@xyLineThicknessF 		= 3.0

res@trYReverse			= True

do ivar=0,dimsizes(var2D)-1

res@gsnLeftString		= var2D(ivar)

res@xyLineColors		= (/"red","blue"/)
res@xyDashPattern		= 0

res@gsnYRefLine			= (/0/)

;plot(0)			 	= gsn_csm_y(wks,vars_LatLon_mm(ivar,:,{lonS:lonE}),res)

tmp1				= vars_LatLon_mm(ivar,:,{lonS:lonE})
plot(0)			 	= gsn_csm_xy(wks,tmp1&lon,tmp1,res)


;---------------------------------------------------------------
; Draw Lon-Value Plots: differences
;---------------------------------------------------------------

res@xyDashPattern		= 0		; solid line

;plot(1)		 		= gsn_csm_y(wks,vars_LatLon_diff_mm(ivar,{lonS:lonE}),res)
tmp2				= vars_LatLon_diff_mm(ivar,{lonS:lonE})
plot(1)		 		= gsn_csm_xy(wks,tmp2&lon,tmp2,res)


; add legends

legends                         = (/"DEF","NEW"/)

lgres                           = True
lgres@lgLineThicknessF 		= res@xyLineThicknessF
lgres@lgLineColors	 	= res@xyLineColors
lgres@lgLabelFontHeightF        = 0.08
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.1
lgres@lgPerimOn                 = False

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)
lbid2                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres                           = True
amres@amJust			= "TopLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = -0.50

annoid1                         = gsn_add_annotation(plot(0),lbid1,amres)
annoid2                         = gsn_add_annotation(plot(1),lbid2,amres)

gsn_panel(wks,plot,(/1,2/),False)

end do ; do ivar=

delete(res)
delete(plot)
delete(wks)
delete(tmp1)
delete(tmp2)

end if ; LonValuePlot


;---------------------------------------------------------------
; Draw Lat-Value Plots
;---------------------------------------------------------------
if(LatValuePlot)then
;---------------------------------------------------------------
; Draw Lon-Value Plots: separate schemes
;---------------------------------------------------------------
wks = gsn_open_wks(PlotType,"./figure/ocn-LatValue_"+lonS+"_"+lonE)

plot				=new(2,graphic)

res 				= True

res@gsnFrame		 	= False
res@gsnDraw			= False

res@xyLineThicknessF 		= 3.0

res@trYReverse			= True

do ivar=0,dimsizes(var2D)-1

res@gsnLeftString		= var2D(ivar)

res@xyLineColors		= (/"red","blue"/)
res@xyDashPattern		= 0

res@gsnYRefLine			= (/0/)

tmp1				= vars_LatLon_zm(ivar,:,{latS:latE})
plot(0)			 	= gsn_csm_xy(wks,tmp1&lat,tmp1,res)

;---------------------------------------------------------------
; Draw Lon-Value Plots: differences
;---------------------------------------------------------------

res@xyDashPattern		= 0		; solid line

tmp2				= vars_LatLon_diff_zm(ivar,{latS:latE})
plot(1)		 		= gsn_csm_xy(wks,tmp2&lat,tmp2,res)

; add legends

legends                         = (/"DEF","NEW"/)

lgres                           = True
lgres@lgLineThicknessF 		= res@xyLineThicknessF
lgres@lgLineColors	 	= res@xyLineColors
lgres@lgLabelFontHeightF        = 0.08
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.1
lgres@lgPerimOn                 = False

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)
lbid2                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres                           = True
amres@amJust			= "TopLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = -0.50

annoid1                         = gsn_add_annotation(plot(0),lbid1,amres)
annoid2                         = gsn_add_annotation(plot(1),lbid2,amres)

gsn_panel(wks,plot,(/1,2/),False)

end do ; do ivar=

delete(res)
delete(plot)
delete(wks)

end if ; LatValuePlot



if(LatLonPlot)then

;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------
do ilev = 0,dimsizes(lev)-1

wks = gsn_open_wks(PlotType,"./figure/ocn-lat-lon")

plot 				= new((/dimsizes(case),dimsizes(var3D)/),graphic)
plotdiff 			= new((/dimsizes(var3D)/),graphic)


res 				= True

res@gsnFrame			= False
res@gsnDraw			= False

res@gsnSpreadColors 		= True

res@cnFillOn 			= True
res@cnLinesOn 			= True
res@cnLineLabelsOn 		= False

res@lbOrientation 		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= -0.01

res@lbLabelStride    		= 2                        ; every other label
;res@tmXBTickSpacingF 		= 60.                      ; every 60 degrees 
;res@tmXBTickSpacingF 		= 15.


res@tmXBLabelFontHeightF  	= 0.014
res@tmYLLabelFontHeightF  	= 0.014

res@mpCenterLonF		= 180.
res@mpMinLatF			= latS
res@mpMaxLatF			= latE
res@mpMinLonF			= lonS
res@mpMaxLonF			= lonE

res@gsnAddCyclic		= False

res@gsnScalarContour     	= True               ; contours desired

res@vcRefMagnitudeF		= 1.0             ; define vector ref mag
res@vcRefLengthF		= 0.045           ; define length of vec ref
res@vcRefAnnoOrthogonalPosF	= -1.0            ; move ref vector
res@vcRefAnnoArrowLineColor	= "black"         ; change ref vector color
res@vcRefAnnoArrowUseVecColor	= False           ; don't use vec color for ref

;res@vcGlyphStyle		= "CurlyVector"     ; turn on curly vectors
res@vcLineArrowColor		= "white"           ; change vector color
res@vcLineArrowThicknessF	= 4.0               ; change vector thickness
res@vcVectorDrawOrder		= "PostDraw"        ; draw vectors last
res@vcMinDistanceF		= 0.025

res@vcMinMagnitudeF		= 1.0		; NOTE: this cannot work when "CurlyVector" is turned on.


string_cli 			= (/"DEF","NEW"/)

resP 				= True
;resP@gsnPanelRowSpec 		= True


;---------------------------------------------------------
; plot figures of separate schemes
;---------------------------------------------------------

if(latE.lt.0)then
res@gsnRightString 		= abs(latS)+"S-"+abs(latE)+"S"
else
res@gsnRightString 		= abs(latS)+"S-"+abs(latE)+"N"
end if


do ivar = 0,dimsizes(var3D)-2
	res@gsnCenterString 		= var3D(ivar)

	maxval				= max(vars_LevLatLon(ivar,:,{lev(ilev)},{latS:latE},{lonS:lonE}))
	minval				= min(vars_LevLatLon(ivar,:,{lev(ilev)},{latS:latE},{lonS:lonE}))
	spcval				= (maxval-minval)/11

	if(var3D(ivar).eq."UVEL".or.var3D(ivar).eq."VVEL".or.var3D(ivar).eq."WVEL")then
		tmp				= (/abs(maxval),abs(minval)/)
		exmval				= min(tmp)
		data_levels			= fspan(-1.*exmval,exmval,11)
	else
		data_levels			= fspan(minval,maxval,11)
	end if
	
	res@cnLevelSelectionMode 	= "ExplicitLevels"
	res@cnLevels             	= data_levels

	res@cnLineLabelFormat    	= "0@;*.3f"                      ; Only one value after decimal point
	if(abs(spcval).lt.1.e-4)then
		res@lbLabelStrings       	= sprintf("%5.1g",data_levels)   ; Format the labelbar labels
	else
		res@lbLabelStrings       	= sprintf("%5.2f",data_levels)   ; Format the labelbar labels
	end if


	do icase = 0,1
		res@gsnLeftString 		= string_cli(icase)

		plot(icase,ivar)	= gsn_csm_vector_scalar_map_ce(wks,vars_LatLon(0,icase,{latS:latE},{lonS:lonE}),vars_LatLon(1,icase,{latS:latE},{lonS:lonE}),vars_LevLatLon(ivar,icase,{lev(ilev)},{latS:latE},{lonS:lonE}),res)


	end do ; do icase=
;	gsn_panel(wks,plot(:,ivar),(/2,1/),resP)
end do ; do ivar=

delete(ivar)
delete(icase)


;---------------------------------------------------------
; plot difference figures
;---------------------------------------------------------

string_diff 			= (/"NEW-DEF"/)
res@gsnLeftString 		= string_diff

res@vcMinMagnitudeF		= 0.20		; NOTE: this cannot work when "CurlyVector" is turned on.
res@vcRefMagnitudeF		= 0.5             ; define vector ref mag

do ivar = 0,dimsizes(var3D)-1
	maxval				= max(vars_LevLatLon_diff(ivar,{lev(ilev)},{latS:latE},{lonS:lonE}))
	minval				= min(vars_LevLatLon_diff(ivar,{lev(ilev)},{latS:latE},{lonS:lonE}))
	spcval				= (maxval-minval)/11

	tmp				= (/abs(maxval),abs(minval)/)
	exmval				= max(tmp)
	data_levels			= fspan(-1.*exmval,exmval,11)
	
	res@cnLevelSelectionMode 	= "ExplicitLevels"
	res@cnLevels             	= data_levels

	res@cnLineLabelFormat    	= "0@;*.3f"                      ; Only one value after decimal point
	if(abs(spcval).lt.1.e-2)then
		res@lbLabelStrings       	= sprintf("%5.1g",data_levels)   ; Format the labelbar labels
	else
		res@lbLabelStrings       	= sprintf("%5.2f",data_levels)   ; Format the labelbar labels
	end if

	res@gsnCenterString 		= var3D(ivar)


	if(var3D(ivar).ne."WVEL")then
		plotdiff(ivar)	= gsn_csm_vector_scalar_map_ce(wks,vars_LatLon_diff(0,{latS:latE},{lonS:lonE}),vars_LatLon_diff(1,{latS:latE},{lonS:lonE}),vars_LevLatLon_diff(ivar,{lev(ilev)},{latS:latE},{lonS:lonE}),res)
	end if

end do ; do ivar=

gsn_panel(wks,plotdiff,(/2,2/),resP)

delete(res)
delete(plot)
delete(wks)

end do ; do ilev=

end if ; LatLonPlot

end do ; do is=
end do ; do icomp=

 
                               
end
