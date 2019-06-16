load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;========================================================
; Read data from my output file: terms.nc
;========================================================
season			 	= "ANN"


ncin			= addfile("terms_"+season+".nc","r")
finalvar_Wm2		= ncin->finalvar_Wm2
finalvardiff_Wm2	= ncin->finalvardiff_Wm2

printVarSummary(finalvar_Wm2)
printVarSummary(finalvardiff_Wm2)


;========================================================
; Prepare plotting data
;========================================================
latS				= -20.
latE				= -5.

lonS				= 160.
lonE				= 280.
;lonS				= 170.
;lonE				= 210.



terms				= (/"residual","heat flux","u_dTdx","v_dTdy","entrainment","sumall"/)

LonValuePlot			= True

LatLonPlot			= False

PlotType			= "x11"

;========================================================
; Heat budget output in one specific domain
; to compare with the results from Liu et al. (2009)
; which have observation data 
;========================================================

lat				= finalvar_Wm2&lat
lat!0				= "lat"
lat&lat				= finalvar_Wm2&lat
printVarSummary(lat)

rad    				= 4.0*atan(1.0)/180.0
clat   				= cos(lat*rad)
copy_VarCoords(lat,clat)

finalvar_Wm2_aave 		= wgt_areaave_Wrap(finalvar_Wm2(:,:,{latS:latE},{lonS:lonE}), clat({latS:latE}), 1.0, 1)
printVarSummary(finalvar_Wm2_aave)

header 				= (/"residual		netFlux		u_dTdx		v_dTdy		entrain"/)
hlist				= [/header/]
alist				= [/finalvar_Wm2_aave(0,:),finalvar_Wm2_aave(1,:),\
				    finalvar_Wm2_aave(2,:),finalvar_Wm2_aave(3,:),\
				    finalvar_Wm2_aave(4,:)/]

fname_t				= "Table_OcnHeatBudget.txt"

write_table(fname_t,"w",hlist,"%s")
write_table(fname_t,"a",alist,"%10.3f%10.3f%10.3f%10.3f%10.3f")



;========================================================
; Get meridional mean over (5S-10S) zonal band 
;========================================================
finalvar_Wm2_mm		= dim_avg_n_Wrap(finalvar_Wm2(:,:,{latS:latE},:),2)
finalvardiff_Wm2_mm	= dim_avg_n_Wrap(finalvardiff_Wm2(:,{latS:latE},:),1)

tmp1 			= finalvar_Wm2(terms|:,case|:,lon|:,lat|:)
copy_VarCoords_1(tmp1,finalvar_Wm2_mm)

tmp2 			= finalvardiff_Wm2(terms|:,lon|:,lat|:)
copy_VarCoords_1(tmp2,finalvardiff_Wm2_mm)

finalvar_Wm2_mm&lon@units		= "degrees_east"
finalvardiff_Wm2_mm&lon@units		= "degrees_east"

print("")
printVarSummary(finalvar_Wm2_mm)
printMinMax(finalvar_Wm2_mm,True)
print("")
printVarSummary(finalvardiff_Wm2_mm)
printMinMax(finalvardiff_Wm2_mm,True)
print("")


;========================================================
;========================================================
; Start plot figures here
;========================================================
;========================================================



if(LonValuePlot)then
;---------------------------------------------------------------
; Draw Lon-Value Plots: separate schemes
;---------------------------------------------------------------
if(latE.lt.0)then
wks = gsn_open_wks(PlotType,"./figure/ocn-heat-budget-LonValue-"+abs(latS)+"S-"+abs(latE)+"S"+"_"+season)
else
wks = gsn_open_wks(PlotType,"./figure/ocn-heat-budget-LonValue-"+abs(latS)+"S-"+abs(latE)+"N"+"_"+season)
end if


plot				=new(2,graphic)

resL 				= True

resL@vpWidthF			= 0.7
resL@vpHeightF			= 0.4

resL@gsnFrame		 	= False
resL@gsnDraw			= False

resL@xyLineThicknessF 		= 3.0

;resL@tmXBMode                   = "Manual"    
;resL@tmXBTickStartF             = lonS+1
;resL@tmXBTickEndF               = lonE-1
resL@tmXBTickSpacingF           = 15.

resL@gsnYRefLine       	= 0.0
resL@gsnYRefLineColor		= "grey" 
resL@gsnYRefLineThicknessF	= 3.0

resL@tiYAxisString		= "Tendency (W/m~S~2~N~)"
;resL@tiXAxisString		= "longitude"


;resL@xyDashPatterns		= (/0,0,0,0,0/)
resL@xyLineColors		= (/"black","red","blue","green","orange","yellow"/)

;plot(0)			 	= gsn_csm_y(wks,finalvar_Wm2_mm(:,0,{lonS:lonE}),resL)
;
;resL@xyDashPattern		= 2
;plot1			 	= gsn_csm_y(wks,finalvar_Wm2_mm(:,1,{lonS:lonE}),resL)
;overlay(plot(0),plot1)

;---------------------------------------------------------------
; Draw Lon-Value Plots: differences
;---------------------------------------------------------------

resL@trYMinF			= round(min(finalvardiff_Wm2_mm(:,{lonS:lonE})),3)-5.
resL@trYMaxF			= round(max(finalvardiff_Wm2_mm(:,{lonS:lonE})),3)+10.

;resL@xyDashPatterns		= (/0,0,0,0,0/)

tmp				= finalvardiff_Wm2_mm(1:,{lonS:lonE})

plot(1)		 		= gsn_csm_xy(wks,tmp&lon,tmp,resL)

;----------------------------------------
; add legends
;----------------------------------------

legends                         = terms(1:)

lgres                           = True
lgres@lgLineColors            	= resL@xyLineColors
lgres@lgLineThicknessF		= resL@xyLineThicknessF
;lgres@lgDashIndexes		= resL@xyDashPatterns
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

;annoid1                         = gsn_add_annotation(plot(0),lbid1,amres)
annoid2                         = gsn_add_annotation(plot(1),lbid2,amres)

gsn_panel(wks,plot(1),(/1,1/),False)


;gsn_panel(wks,plot,(/1,2/),False)

delete(resL)
delete(plot)
delete(wks)

end if ; LonValuePlot


;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------

if(LatLonPlot)then

wks = gsn_open_wks(PlotType,"./figure/ocn-heat-budget-LatLon")

plot = new(5,graphic)

lat_min 			= latS
lat_max 			= latE
lon_min 			= lonS
lon_max 			= lonE

res 				= True
res@gsnFrame 			= False
res@gsnDraw 			= False
res@gsnSpreadColors 		= True

res@cnFillOn 			= True
res@cnLinesOn 			= False
res@cnLineLabelsOn 		= False

res@mpLandFillColor 		= "gray"
res@mpMinLonF 			= lon_min
res@mpMaxLonF 			= lon_max
res@mpMinLatF 			= lat_min
res@mpMaxLatF 			= lat_max

res@mpCenterLonF 		= 180.

res@lbOrientation 		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= -0.01

res@lbLabelStride    		= 2                        ; every other label
res@tmXBTickSpacingF 		= 60.                      ; every 60 degrees 

;res@mpGridAndLimbOn        	= True               ; turn on grid lines
;res@mpGridLineDashPattern 	= 5
;res@mpGridLatSpacingF 		= 10.
;res@mpGridLonSpacingF 		= 20.
;res@mpGridMaskMode		= "MaskNotOcean"

res@tmXBLabelFontHeightF  	= 0.014
res@tmYLLabelFontHeightF  	= 0.014

res@gsnLeftStringFontHeightF	= 0.030
res@gsnRightStringFontHeightF	= 0.030
res@gsnCenterStringFontHeightF	= 0.030


string_cli 			= (/"DEF","NEW"/)

resP 				= True
resP@gsnPanelRowSpec 		= True


;---------------------------------------------------------
; plot figures with units of Wm2
;---------------------------------------------------------

res@cnLevelSelectionMode 	= "ManualLevels" ; manual levels
res@cnMinLevelValF       	= -100.
res@cnMaxLevelValF       	= 100.
res@cnLevelSpacingF      	= 10.

res@gsnRightString 		= "W/m~S~2~N~"

do icli = 0,1
res@gsnLeftString 		= string_cli(icli)
do ifin = 0,4
res@gsnCenterString 		= terms(ifin)
	plot(ifin) 		= gsn_csm_contour_map_ce(wks,finalvar_Wm2(ifin,icli,:,:),res) 

end do ; do ifin=
gsn_panel(wks,plot,(/1,2,2/),resP)
end do ; do icli=

; diff
res@cnLevelSelectionMode 	= "ManualLevels" ; manual levels
res@cnMinLevelValF       	= -50.
res@cnMaxLevelValF       	= 50.
res@cnLevelSpacingF      	= 5.

string_diff 			= (/"NEW-DEF"/)
res@gsnLeftString 		= string_diff

do ifin = 0,4
res@gsnCenterString 		= terms(ifin)
	plot(ifin) 		= gsn_csm_contour_map_ce(wks,finalvardiff_Wm2(ifin,:,:),res) 

end do ; do ifin=


;---------------------------------------------------------------
; add one box over eastern pacific that we care more about.
;---------------------------------------------------------------
dum1 				= new(4,graphic)

resp 				= True
resp@gsLineColor 		= "blue"
resp@gsLineThicknessF 		= 3.0

ypts 				= (/30.0, 30.0, -30.0, -30.0, 30.0/)
xpts 				= (/180.0, 220.0, 220.0, 180.0, 180.0/)

; draw each ine separately. Each line must contain two points.
;do i = 0,3
;	dum1(i) 		= gsn_add_polyline(wks,plot(0),xpts(i:i+1),ypts(i:i+1),resp)
;end do


gsn_panel(wks,plot,(/1,2,2/),resP)

delete(res)
delete(plot)
delete(wks)

end if ; LatLonPlot
                               
end
