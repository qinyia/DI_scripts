load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************

TEMP_levs		= fspan(0,30,31)
TEMP_dlevs		= fspan(-1.5,1.5,31)

UVEL_levs		= fspan(-30,30,31)
UVEL_dlevs		= fspan(-15,15,31)


VVEL_levs		= fspan(-7.5,7.5,31)
VVEL_dlevs		= fspan(-3,3,31)

season			=(/"ANN"/)

var3D			=(/"UVEL","VVEL","TEMP"/)
var3D_levs		= (/UVEL_levs,VVEL_levs,TEMP_levs/)
var3D_dlevs		= (/UVEL_dlevs,VVEL_dlevs,TEMP_dlevs/)
var3D_units		= (/"cm/s","cm/s","degC"/)


var2D			=(/"HMXL","SHF"/)
compare_type 		= (/"cmip-obs"/)

Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do is = 0,dimsizes(season)-1

ncdf = addfile("ocn-cross-section-w-obs-vars_3D_"+season(is)+".nc" ,"r")

dSizes			= getfiledimsizes(ncdf)

vars_3D			= new((/dimsizes(var3D),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_3D!0		= "vars"
vars_3D&vars		= ispan(0,dimsizes(var3D)-1,1)

do ivar=0,dimsizes(var3D)-1
	vars_3D(ivar,:,:,:,:)		= ncdf->$var3D(ivar)$

; change units from m/s to cm/s
	if(var3D(ivar).eq."UVEL".or.var3D(ivar).eq."VVEL")then
		vars_3D(ivar,:,:,:,:) 	= vars_3D(ivar,:,:,:,:)*1e2
	end if

end do ; do ivar=

printVarSummary(vars_3D)
;========================================================
; Get the difference between two cases
;========================================================
vars_3D_diff			= vars_3D
vars_3D_diff(:,0,:,:,:)		= vars_3D(:,1,:,:,:)-vars_3D(:,0,:,:,:)
vars_3D_diff(:,1,:,:,:)		= vars_3D(:,2,:,:,:)-vars_3D(:,0,:,:,:)
vars_3D_diff(:,2,:,:,:)		= vars_3D(:,2,:,:,:)-vars_3D(:,1,:,:,:)

printVarSummary(vars_3D_diff)

; get UVEL and VVEL to plot vector.
U_3D				= vars_3D(0,:,:,:,:)
V_3D				= vars_3D(1,:,:,:,:)
T_3D				= vars_3D(2,:,:,:,:)

U_3D_diff			= vars_3D_diff(0,:,:,:,:)
V_3D_diff			= vars_3D_diff(1,:,:,:,:)
T_3D_diff			= vars_3D_diff(2,:,:,:,:)

printVarSummary(U_3D)
printVarSummary(V_3D_diff)

;========================================================
; Prepare plot settings
;========================================================

latS			= -30.
latE			= 30.

lonS			= 160.
lonE			= 280.

levS			= 50.

lev			= (/6,10,30,50,100,200,300/)
;lev			= (/6/)


;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------
do ilev = 0,dimsizes(lev)-1

wks = gsn_open_wks("pdf","./figure/ocn-current-vector-at-"+lev(ilev)+"m_revised")


plot 				= new((/dimsizes(var3D)/),graphic)
plotdiff 			= new((/dimsizes(var3D)/),graphic)


res 				= True

res@gsnFrame			= False
res@gsnDraw			= False

;res@gsnSpreadColors 		= True

res@gsnStringFontHeightF	= 0.032
res@lbLabelFontHeightF		= 0.032
res@tmXBLabelFontHeightF	= 0.032
res@tmYLLabelFontHeightF	= 0.032

res@mpCenterLonF		= 180.
res@mpMinLatF			= latS
res@mpMaxLatF			= latE
res@mpMinLonF			= lonS
res@mpMaxLonF			= lonE

res@gsnAddCyclic		= False

cmap1				= read_colormap_file("BlAqGrYeOrReVi200")
cmap2				= read_colormap_file("BlueRed")

res@lbLabelStride		= 2
res@vcRefMagnitudeF	 	= 20.0
res@vcRefLengthF		= 0.05
res@vcMinDistanceF	 	= 0.02
res@vcGlyphStyle		= "CurlyVector"
res@vcLevelPalette		= cmap1(6:193,:)

res@vcLevelSelectionMode	= "ExplicitLevels"
res@vcLevels			= TEMP_levs

res@vcLineArrowThicknessF	= 1.5

res@vcLineArrowHeadMinSizeF	= 0.01
res@vcLineArrowHeadMaxSizeF	= 0.01

string_cli 			= (/"OBS","DEF","NEW"/)

resP 				= True
resP@gsnPanelMainString		= lev(ilev)+" m"

resP@amJust                     = "TopRight"
resP@gsnPanelFigureStringsFontHeightF	= 0.02
resP@gsnPanelFigureStrings      = (/"a)","b)","c)","d)","e)","f)"/)

;---------------------------------------------------------
; plot figures of separate schemes
;---------------------------------------------------------

if(latE.lt.0)then
	res@gsnRightString 		= abs(latS)+"S-"+abs(latE)+"S"
else
	res@gsnRightString 		= abs(latS)+"S-"+abs(latE)+"N"
end if

do icase = 0,2
	res@gsnLeftString 		= string_cli(icase)
	U_tmp				= U_3D(icase,{lev(ilev)},{latS:latE},{lonS:lonE})
	V_tmp				= V_3D(icase,{lev(ilev)},{latS:latE},{lonS:lonE})
	T_tmp				= T_3D(icase,{lev(ilev)},{latS:latE},{lonS:lonE})
	plot(icase)	= gsn_csm_vector_scalar_map_ce(wks,U_tmp(::3,::3),V_tmp(::3,::3),T_tmp(::3,::3),res) 

end do ; do icase=

;gsn_panel(wks,plot,(/1,3/),resP)

delete(icase)


;---------------------------------------------------------
; plot difference figures
;---------------------------------------------------------

res@vcRefMagnitudeF	 	= 7.0
res@vcRefLengthF		= 0.05

res@vcLevelSelectionMode	= "ExplicitLevels"
res@vcLevels			= TEMP_dlevs

delete(res@vcLevelPalette)
res@vcLevelPalette		= cmap2


string_diff 			= (/"DEF-OBS","NEW-OBS","NEW-DEF"/)

do icase = 0,dimsizes(string_diff)-1
	res@gsnLeftString 		= string_diff(icase)

	U_tmp				= U_3D_diff(icase,{lev(ilev)},{latS:latE},{lonS:lonE})
	V_tmp				= V_3D_diff(icase,{lev(ilev)},{latS:latE},{lonS:lonE})
	T_tmp				= T_3D_diff(icase,{lev(ilev)},{latS:latE},{lonS:lonE})
	plotdiff(icase)	= gsn_csm_vector_scalar_map_ce(wks,U_tmp(::3,::3),V_tmp(::3,::3),T_tmp(::3,::3),res) 
	end do ; do icase=

;gsn_panel(wks,plotdiff,(/1,3/),resP)

plotall		= (/plot(0),plotdiff(0),plot(1),plotdiff(1),plot(2),plotdiff(2)/)
gsn_panel(wks,ndtooned(plotall),(/3,2/),resP)
delete(plotall)
plotall		= (/plot,plotdiff/)
gsn_panel(wks,ndtooned(plotall),(/2,3/),resP)
delete(plotall)


delete(res)
delete(plot)
delete(wks)

end do ; do ilev=

end do ; do is=

 
                               
end
