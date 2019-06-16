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

var2D			=(/"U925","V925","SST","PRECT","PSL"/)

compare_type 		= (/"cmip-obs"/)
case			= (/"def","new"/)

Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp=0,dimsizes(compare_type)-1
do is=0,dimsizes(season)-1

; for 2D varialbes
ncdf = addfile("atm-var2D_"+compare_type(icomp)+"_"+season(is)+"_w_obs.nc","r")

dSizes              = getfiledimsizes(ncdf)  

vars_LatLon         = new((/dimsizes(var2D),dSizes(0),dSizes(1),dSizes(2)/),float)

do ivar=0,dimsizes(var2D)-1
        vars_LatLon(ivar,:,:,:)           =ncdf->$var2D(ivar)$ 
end do

vars_LatLon!0		= "vars"
vars_LatLon&lat@units	= "degree_north"
vars_LatLon&lon@units	= "degree_east"

printVarSummary(vars_LatLon)
printMinMax(vars_LatLon,True)

;========================================================
; Prepare plot settings
;========================================================

latS			= -30.
latE			= 30.

lonS			= 160.
lonE			= 280.

levS			= 50.

WindMagPlot		= False
LonValuePlot		= False
LatValuePlot		= False

LatLonPlot		= True

PlotType		= "x11"


;========================================================
; Get the difference between two cases
;========================================================

vars_LatLon_diff			= vars_LatLon
vars_LatLon_diff(:,0,:,:)		= vars_LatLon(:,1,:,:)-vars_LatLon(:,0,:,:)
vars_LatLon_diff(:,1,:,:)		= vars_LatLon(:,2,:,:)-vars_LatLon(:,0,:,:)
vars_LatLon_diff(:,2,:,:)		= vars_LatLon(:,2,:,:)-vars_LatLon(:,1,:,:)
printVarSummary(vars_LatLon_diff)

;========================================================
; Get the wind magnitude
;========================================================
uspeed				= vars_LatLon(0,:,:,:)
vspeed				= vars_LatLon(1,:,:,:)
windmag				= sqrt(uspeed^2+vspeed^2)
copy_VarCoords(uspeed,windmag)
printVarSummary(windmag)

windmag_diff			= windmag
windmag_diff(0,:,:)			= windmag(1,:,:)-windmag(0,:,:)
windmag_diff(1,:,:)			= windmag(2,:,:)-windmag(0,:,:)
windmag_diff(2,:,:)			= windmag(2,:,:)-windmag(1,:,:)
printMinMax(windmag_diff,True)

;========================================================
; Get the meridional mean over specific latitude ranges--> Lon-value plots
;========================================================

vars_LatLon_mm			= dim_avg_n_Wrap(vars_LatLon(:,:,{latS:latE},:),2)
printVarSummary(vars_LatLon_mm)

vars_LatLon_diff_mm		= dim_avg_n_Wrap(vars_LatLon_diff(:,:,{latS:latE},:),2)
printVarSummary(vars_LatLon_diff_mm)
printMinMax(vars_LatLon_diff_mm,True)

;========================================================
; Get the zonal mean over specific longitude ranges--> Lat-Value plots
;========================================================
vars_LatLon_zm			= dim_avg_n_Wrap(vars_LatLon(:,:,:,{lonS:lonE}),3)
printVarSummary(vars_LatLon_zm)

vars_LatLon_diff_zm		= dim_avg_n_Wrap(vars_LatLon_diff(:,:,:,{lonS:lonE}),3)
printVarSummary(vars_LatLon_diff_zm)
printMinMax(vars_LatLon_diff_zm,True)

;--------------------------------------------------------------------------------------------------------
; Plotting starts here...
;--------------------------------------------------------------------------------------------------------

if(WindMagPlot)then

;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------
wks = gsn_open_wks(PlotType,"./figure/atm-lat-lon-windmag-overlay-SST")

gsn_define_colormap(wks,"BlueDarkRed18")

plot 				= new(3,graphic)
plotdiff			= new(3,graphic)


res 				= True

res@gsnFrame			= False
res@gsnDraw			= False

res@gsnSpreadColors 		= True

res@cnFillOn 			= True
res@cnLinesOn 			= False
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

res@gsnRightString		= "SST (shaded) & wind speed (contour)"

; overlay resources
resov			 	= True
resov@cnFillOn			= False
resov@cnLinesOn			= True
resov@cnLineLabelsOn		= False
resov@gsnFrame			= False
resov@gsnDraw			= False

resov@cnLineThicknessF		= 3.0

resov@cnInfoLabelOn   		= False

resov@gsnContourNegLineDashPattern = 1
resov@gsnContourPosLineDashPattern = 0
;resov@gsnContourZeroLineThicknessF = 0.0



string_cli 			= (/"OBS","DEF","NEW"/)

resP 				= True

;---------------------------------------------------------
; plot figures of separate schemes
;---------------------------------------------------------
res@cnLevelSelectionMode	= "ManualLevels"
res@cnMinLevelValF		= 18.
res@cnMaxLevelValF		= 30.
res@cnLevelSpacingF		= 1.

resov@cnLevelSelectionMode	= "ManualLevels"
resov@cnMinLevelValF		= 1.
resov@cnMaxLevelValF		= 11.
resov@cnLevelSpacingF		= 2.



do icase = 0,2
	res@gsnLeftString 		= string_cli(icase)
	SST_tmp				= vars_LatLon(2,icase,{latS:latE},{lonS:lonE})
	windmag_tmp			= windmag(icase,{latS:latE},{lonS:lonE})

	plot(icase)	= gsn_csm_contour_map(wks,SST_tmp,res)
	plot_ov		= gsn_csm_contour(wks,windmag_tmp,resov)
	overlay(plot(icase),plot_ov)

end do ; do icase=

;gsn_panel(wks,plot(0:2),(/3,1/),resP)

delete(SST_tmp)
delete(windmag_tmp)

;---------------------------------------------------------
; plot difference figures
;---------------------------------------------------------
res@cnLevelSelectionMode	= "ManualLevels"
res@cnMinLevelValF		= -5.
res@cnMaxLevelValF		= 5.
res@cnLevelSpacingF		= 1.

resov@cnLevelSelectionMode	= "ManualLevels"
resov@cnMinLevelValF		= -5.
resov@cnMaxLevelValF		= 5.
resov@cnLevelSpacingF		= 1.


string_diff 			= (/"NEW-OBS","DEF-OBS","NEW-DEF"/)

do icase=0,2
res@gsnLeftString 		= string_diff(icase)

SST_tmp				= vars_LatLon_diff(2,icase,{latS:latE},{lonS:lonE})
windmag_tmp			= windmag_diff(icase,{latS:latE},{lonS:lonE})

plotdiff(icase)				= gsn_csm_contour_map(wks,SST_tmp,res)
plotdiff_ov				= gsn_csm_contour(wks,windmag_tmp,resov)
overlay(plotdiff(icase),plotdiff_ov)

end do

plotall		= new(6,graphic)
plotall(::2)	= plot
plotall(1::2)	= plotdiff

;gsn_panel(wks,plotdiff,(/3,1/),resP)
gsn_panel(wks,plotall,(/3,2/),resP)

delete(res)
delete(plot)
delete(wks)


end if ; WindCurlPlot



;---------------------------------------------------------------
; Draw Lon-Value Plots
;---------------------------------------------------------------
if(LonValuePlot)then

latSs			= (/-15,-5,5/)
latEs			= (/-5,5,15/)


wks = gsn_open_wks(PlotType,"./figure/atm-LonValue")

plot				= new((/dimsizes(latSs),dimsizes(var2D),2/),graphic)

do ilat=0,dimsizes(latSs)-1

;========================================================
; Get the meridional mean over specific latitude ranges--> Lon-value plots
;========================================================

vars_LatLon_mm_S			= dim_avg_n_Wrap(vars_LatLon(:,:,{latSs(ilat):latEs(ilat)},:),2)

vars_LatLon_diff_mm_S			= dim_avg_n_Wrap(vars_LatLon_diff(:,:,{latSs(ilat):latEs(ilat)},:),2)

;---------------------------------------------------------------
; Draw Lon-Value Plots: separate schemes
;---------------------------------------------------------------

;plot				=new(2,graphic)

res 				= True

res@gsnFrame		 	= False
res@gsnDraw			= False

res@xyLineThicknessF 		= 3.0

;res@trYReverse			= True

res@vpWidthF			= 0.75
res@vpHeightF			= 0.40

do ivar=0,dimsizes(var2D)-1

res@gsnLeftString		= var2D(ivar)

res@xyLineColors		= (/"grey","red","blue"/)
res@xyDashPattern		= 0

;res@gsnYRefLine			= (/0/)

res@gsnRightString		= latSs(ilat)+"_"+latEs(ilat)


tmp1						= vars_LatLon_mm_S(ivar,:,{lonS:lonE})
plot(ilat,ivar,0)			 	= gsn_csm_xy(wks,tmp1&lon,tmp1,res)


;---------------------------------------------------------------
; Draw Lon-Value Plots: differences
;---------------------------------------------------------------

res@xyDashPattern		= 0		; solid line

tmp2						= vars_LatLon_diff_mm_S(ivar,:,{lonS:lonE})
plot(ilat,ivar,1)		 		= gsn_csm_xy(wks,tmp2&lon,tmp2,res)


; add legends

legends                         = (/"OBS","DEF","NEW"/)
legendsDiff                         = (/"DEF-OBS","NEW-OBS","NEW-DEF"/)

lgres                           = True
lgres@lgLineThicknessF 		= res@xyLineThicknessF
lgres@lgLineColors	 	= res@xyLineColors
lgres@lgDashIndexes		= (/0,0,0/)
lgres@lgLabelFontHeightF        = 0.08
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.1
lgres@lgPerimOn                 = False

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)
lbid2                           = gsn_create_legend(wks,dimsizes(legendsDiff),legendsDiff,lgres)

amres                           = True
amres@amJust			= "BottomLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = 0.50

if(ivar.eq.0)then ; only add legend for the first variable plot
annoid1                         = gsn_add_annotation(plot(ilat,ivar,0),lbid1,amres)
annoid2                         = gsn_add_annotation(plot(ilat,ivar,1),lbid2,amres)
end if


end do ; do ivar=

delete(res)
delete(tmp1)
delete(tmp2)

end do ; do ilat=

gsn_panel(wks,ndtooned(plot(:,:,0)),(/dimsizes(latSs),dimsizes(var2D)/),False)

end if ; LonValuePlot


;---------------------------------------------------------------
; Draw Lat-Value Plots
;---------------------------------------------------------------
if(LatValuePlot)then

lonSs			= (/140,180,240/)
lonEs			= (/180,230,280/)

wks = gsn_open_wks(PlotType,"./figure/atm-LatValue")

plot	= new((/dimsizes(lonSs),dimsizes(var2D),2/),graphic)

do ilon=0,dimsizes(lonSs)-1

;========================================================
; Get the zonal mean over specific longitude ranges--> Lat-Value plots
;========================================================
vars_LatLon_zm_S			= dim_avg_n_Wrap(vars_LatLon(:,:,:,{lonSs(ilon):lonEs(ilon)}),3)

vars_LatLon_diff_zm_S		= dim_avg_n_Wrap(vars_LatLon_diff(:,:,:,{lonSs(ilon):lonEs(ilon)}),3)

;---------------------------------------------------------------
; Draw Lon-Value Plots: separate schemes
;---------------------------------------------------------------

res 				= True

res@gsnFrame		 	= False
res@gsnDraw			= False

res@xyLineThicknessF 		= 3.0

res@trYReverse			= True

res@gsnRightString		= lonSs(ilon)+"_"+lonEs(ilon)

do ivar=0,dimsizes(var2D)-1

res@gsnLeftString		= var2D(ivar)

res@xyLineColors		= (/"grey","red","blue"/)
res@xyDashPattern		= 0


tmp1				= vars_LatLon_zm_S(ivar,:,{latS:latE})
plot(ilon,ivar,0)			 	= gsn_csm_xy(wks,tmp1&lat,tmp1,res)

;---------------------------------------------------------------
; Draw Lon-Value Plots: differences
;---------------------------------------------------------------

res@xyDashPattern		= 0		; solid line

tmp2				= vars_LatLon_diff_zm_S(ivar,:,{latS:latE})
plot(ilon,ivar,1)		 		= gsn_csm_xy(wks,tmp2&lat,tmp2,res)

; add legends

legends                         = (/"OBS","DEF","NEW"/)
legendsDiff			= (/"DEF-OBS","NEW-OBS","DEF-NEW"/)

lgres                           = True
lgres@lgLineThicknessF 		= res@xyLineThicknessF
lgres@lgLineColors	 	= res@xyLineColors
lgres@lgDashIndexes		= (/0,0,0/)
lgres@lgLabelFontHeightF        = 0.08
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.1
lgres@lgPerimOn                 = False

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)
lbid2                           = gsn_create_legend(wks,dimsizes(legendsDiff),legendsDiff,lgres)

amres                           = True
amres@amJust			= "TopLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = -0.50

if(ivar.eq.0)then
annoid1                         = gsn_add_annotation(plot(ilon,ivar,0),lbid1,amres)
annoid2                         = gsn_add_annotation(plot(ilon,ivar,1),lbid2,amres)
end if

end do ; do ivar=

delete(res)
delete(tmp1)
delete(tmp2)

end do ; do ilon=

gsn_panel(wks,ndtooned(plot(:,:,0)),(/dimsizes(lonSs),dimsizes(var2D)/),False)

end if ; LatValuePlot


if(LatLonPlot)then

;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------

wks = gsn_open_wks(PlotType,"./figure/atm-lat-lon-windvector-PSL-overlay-SST")

gsn_define_colormap(wks,"BlueDarkRed18")
;gsn_define_colormap(wks,"CBR_coldhot")
;gsn_define_colormap(wks,"MPL_coolwarm")


plot 				= new(3,graphic)
plotdiff			= new(3,graphic)

vcres				= True
vcres@gsnDraw			= False
vcres@gsnFrame			= False
vcres@gsnStringFontHeightF	= 0.035
vcres@tmXBLabelFontHeightF  	= 0.025
vcres@tmYLLabelFontHeightF  	= 0.025

cnres				= vcres
stres				= vcres
mpres				= vcres

;---Set up some vector resources.
vcres@vcGlyphStyle		= "LineArrow"

;---Vector lengths and color
vcres@vcLineArrowThicknessF    =  1.5

vcres@vcLineArrowColor		= "brown"

vcres@vcRefMagnitudeF		= 10.
vcres@vcRefLengthF		= 0.045

vcres@vcRefAnnoString2On	= False
vcres@vcRefAnnoOrthogonalPosF	= -0.535

;---Make sure vectors are drawn in "predraw" phase.
vcres@vcVectorDrawOrder        = "Predraw"

stres@cnFillOn			= True
stres@cnLinesOn			= False
stres@cnLineLabelsOn		= False

stres@cnLevelSelectionMode = "ManualLevels"
stres@cnMinLevelValF       = 20.
stres@cnMaxLevelValF       = 32.
stres@cnLevelSpacingF      = 1.

stres@lbLabelStride		= 2.
stres@pmLabelBarOrthogonalPosF	= 0.15

stres@cnFillDrawOrder		= "PreDraw"

;---Turn on contour fill, and turn other things off.
cnres@cnFillOn             = False
cnres@cnLinesOn            = True
cnres@cnInfoLabelOn        = False

;---Define contour levels
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMinLevelValF       = 984.0
cnres@cnMaxLevelValF       = 1040.0
cnres@cnLevelSpacingF      = 2.0

cnres@cnLineThicknessF		= 2.0
;cnres@cnLineColor		= "grey"
cnres@gsnContourNegLineDashPattern = 2

;---Make sure contours are drawn in "predraw" phase.
cnres@cnFillDrawOrder      = "Predraw"

;---Control appearance of map.
;mpres@mpProjection           = "LambertEqualArea"
mpres@mpLabelsOn             = False
mpres@mpPerimOn              = True
mpres@mpGridAndLimbOn        = False
mpres@mpFillOn               = True
mpres@mpOutlineOn            = True
mpres@mpOutlineDrawOrder     = "PostDraw"
;mpres@mpFillDrawOrder        = "Predraw"
mpres@mpFillDrawOrder        = "Postdraw"

;mpres@mpOceanFillColor       = "lightskyblue1"
mpres@mpLandFillColor        = "gray"

;---Zoom in on area that is roughly the United States.
mpres@mpLimitMode            = "LatLon"
mpres@mpMinLatF              = latS
mpres@mpMaxLatF              = latE
mpres@mpMinLonF              = lonS
mpres@mpMaxLonF              = lonE
mpres@mpCenterLonF           = 180.

mpres@mpGridAndLimbDrawOrder = "Predraw"
;mpres@tiMainString           = title
;
; Create, but don't draw, a vector, streamline, contour
; and map plot. Stride the vector and contour data to
; thin it a little.
;

string_cli 			= (/"OBS","DEF","NEW"/)

do icase = 0,2
	vcres@gsnLeftString 		= string_cli(icase)
	SST_tmp				= vars_LatLon(2,icase,{latS:latE},{lonS:lonE})
	U925_tmp			= vars_LatLon(0,icase,{latS:latE},{lonS:lonE})
	V925_tmp			= vars_LatLon(1,icase,{latS:latE},{lonS:lonE})
	PSL_tmp				= vars_LatLon(4,icase,{latS:latE},{lonS:lonE})

	vcid				= gsn_csm_vector(wks,U925_tmp(::2,::2),V925_tmp(::2,::2),vcres)
	stid				= gsn_csm_contour(wks,SST_tmp,stres)	
	cnid				= gsn_csm_contour(wks,PSL_tmp,cnres)
	plot(icase)			 = gsn_csm_map(wks,mpres)

;---Overlay contour, streamline, and vector plots on the map plot.
	overlay(plot(icase),stid)
	overlay(plot(icase),cnid)
	overlay(plot(icase),vcid)
;	maximize_output(wks,True)

;;;;;;;;;	delete([/mpid,cnid,vcid/])
	delete([/SST_tmp,U925_tmp,V925_tmp,PSL_tmp/])
end do ; do icase=

;gsn_panel(wks,plot,(/3,1/),False)

;---------------------------------------------------------
; plot difference figures
;---------------------------------------------------------
stres@cnLevelSelectionMode	= "ManualLevels"
stres@cnMinLevelValF		= -3.
stres@cnMaxLevelValF		= 3.
stres@cnLevelSpacingF		= 0.5

vcres@vcRefMagnitudeF          = 5.

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMinLevelValF       = -15.
cnres@cnMaxLevelValF       = 15.
cnres@cnLevelSpacingF      = 0.5


string_diff 			= (/"DEF-OBS","NEW-OBS","NEW-DEF"/)

do icase=0,2
	vcres@gsnLeftString 		= string_diff(icase)

	SST_tmp				= vars_LatLon_diff(2,icase,{latS:latE},{lonS:lonE})
	U925_tmp			= vars_LatLon_diff(0,icase,{latS:latE},{lonS:lonE})
	V925_tmp			= vars_LatLon_diff(1,icase,{latS:latE},{lonS:lonE})
	PSL_tmp				= vars_LatLon_diff(4,icase,{latS:latE},{lonS:lonE})

	vcid				= gsn_csm_vector(wks,U925_tmp(::2,::2),V925_tmp(::2,::2),vcres)
	stid				= gsn_csm_contour(wks,SST_tmp,stres)	
	cnid				= gsn_csm_contour(wks,PSL_tmp,cnres)
	plotdiff(icase)			 = gsn_csm_map(wks,mpres)

;---Overlay contour, streamline, and vector plots on the map plot.
	overlay(plotdiff(icase),stid)
	overlay(plotdiff(icase),cnid)
	overlay(plotdiff(icase),vcid)
;	maximize_output(wks,True)

;	delete([/mpid,cnid,vcid/])
	delete([/SST_tmp,U925_tmp,V925_tmp,PSL_tmp/])
end do

;gsn_panel(wks,plotdiff,(/3,1/),False)

plotall		= new(6,graphic)
plotall(::2)	= plot
plotall(1::2)	= plotdiff

resP				= True
resP@amJust			= "TopLeft"
resP@gsnPanelFigureStrings	= (/"a)","b)","c)"/)


;gsn_panel(wks,plotdiff,(/2,2/),resP)
;gsn_panel(wks,plotdiff,(/1,3/),resP)
;gsn_panel(wks,plotdiff,(/3,1/),resP)

resP1				= True
resP1@amJust			= "TopLeft"
resP1@gsnPanelFigureStrings	= (/"a)","b)","c)","d)","e)","f)"/)
gsn_panel(wks,plotall,(/3,2/),resP1)

delete(plot)
delete(wks)

end if ; LatLonPlot
exit

end do ; do is=
end do ; do icomp=

 
                               
end
