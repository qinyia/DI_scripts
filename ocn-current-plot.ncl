load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************

;TEMP_levs		= (/-2,0,2,4,6,10,14,18,20,22,24,26,28/)
TEMP_levs		= fspan(0,30,31)
TEMP_dlevs		= fspan(-3,3,31)
TEMP_dlevs_sp		= (/-3,-2.5,-2,-1.5,-1,1,1.5,2,2.5,3/)


;UVEL_levs		= (/-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12/)
UVEL_levs		= fspan(-30,30,31)
UVEL_dlevs		= fspan(-15,15,31)
;UVEL_dlevs_sp		= (/-5.0,-4.0,-3.0,-2.0,-1.0,1.0,2.0,3.0,4.0,5.0/)


;VVEL_levs		= (/-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0,1.0,2.0,3.0,4.0,5.0,6.0/)
VVEL_levs		= fspan(-7.5,7.5,31)
VVEL_dlevs		= fspan(-3,3,31)
VVEL_dlevs_sp		= (/-1.8,-1.5,-1.2,-0.9,-0.6,0.6,0.9,1.2,1.5,1.8/)

VEL_levs		= (/-2,-1.5,-1.0,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,1.0,1.5,2.0/)
WVEL_dlevs		= (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6/)
WVEL_dlevs_sp		= (/-0.6,-0.5,-0.4,-0.3,-0.2,0.2,0.3,0.4,0.5,0.6/)


;season=(/"ANN","DJF","MAM","JJA","SON"/)
season			=(/"ANN"/)

var3D			=(/"TEMP","UVEL","VVEL"/)
var3D_levs		= (/TEMP_levs,UVEL_levs,VVEL_levs/)
var3D_dlevs		= (/TEMP_dlevs,UVEL_dlevs,VVEL_dlevs/)
var3D_units		= (/"degC","cm/s","cm/s"/)


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

;========================================================
; Prepare plot settings
;========================================================

latS			= -50.
latE			= 50.

lonS			= 160.
lonE			= 280.

levS			= 50.

;lev			= (/10,30,50,100,200,300,500,800,1000/)
lev			= (/6/)


;---------------------------------------------------------------
; Draw Lat-Lon Plots
;---------------------------------------------------------------
do ilev = 0,dimsizes(lev)-1

wks = gsn_open_wks("pdf","./figure/ocn-current-at-"+lev(ilev)+"m_revised")

gsn_define_colormap(wks,"BlueWhiteOrangeRed")

plot 				= new((/3,dimsizes(var3D)/),graphic)
plotdiff 			= new((/3,dimsizes(var3D)/),graphic)


res 				= True

res@gsnFrame			= False
res@gsnDraw			= False

res@gsnSpreadColors 		= True

res@cnFillOn 			= True
res@cnLinesOn 			= False
res@cnLineLabelsOn 		= False

res@gsnStringFontHeightF	= 0.032
res@lbLabelFontHeightF		= 0.032
res@tmXBLabelFontHeightF	= 0.032
res@tmYLLabelFontHeightF	= 0.032


;res@lbOrientation 		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= 0.2

res@lbLabelStride    		= 2                        ; every other label
;res@tmXBTickSpacingF 		= 60.                      ; every 60 degrees 
;res@tmXBTickSpacingF 		= 15.


res@mpCenterLonF		= 180.
res@mpMinLatF			= latS
res@mpMaxLatF			= latE
res@mpMinLonF			= lonS
res@mpMaxLonF			= lonE

res@gsnAddCyclic		= False


string_cli 			= (/"OBS","DEF","NEW"/)

resP 				= True
resP@gsnPanelMainString		= lev(ilev)+" m"


;---------------------------------------------------------
; plot figures of separate schemes
;---------------------------------------------------------

if(latE.lt.0)then
res@gsnRightString 		= abs(latS)+"S-"+abs(latE)+"S"
else
res@gsnRightString 		= abs(latS)+"S-"+abs(latE)+"N"
end if


do ivar = 0,dimsizes(var3D)-1
	res@gsnCenterString 		= var3D(ivar)

	res@cnLevelSelectionMode 	= "ExplicitLevels"
	res@cnLevels             	= var3D_levs(ivar,:)

	do icase = 0,2
	print("ivar="+ivar+", icase="+icase)
		res@gsnLeftString 		= string_cli(icase)
		print(lev(ilev))
		plot(icase,ivar)	= gsn_csm_contour_map_ce(wks,vars_3D(ivar,icase,{lev(ilev)},{latS:latE},{lonS:lonE}),res) 

	end do ; do icase=
;	gsn_panel(wks,plot(:,ivar),(/3,1/),resP)
end do ; do ivar=

plotall	= (/plot(:,0),plot(:,1),plot(:,2)/)
gsn_panel(wks,ndtooned(plotall),(/3,3/),resP)
delete(plotall)

delete(ivar)
delete(icase)

delete(res@cnLevels)


;---------------------------------------------------------
; plot difference figures
;---------------------------------------------------------

string_diff 			= (/"DEF-OBS","NEW-OBS","NEW-DEF"/)

do ivar = 0,dimsizes(var3D)-1

	res@cnLevelSelectionMode 	= "ExplicitLevels"
	res@cnLevels             	= var3D_dlevs(ivar,:)

	res@gsnCenterString 		= var3D(ivar)

	do icase = 0,dimsizes(string_diff)-1
		res@gsnLeftString 		= string_diff(icase)
		plotdiff(icase,ivar)	= gsn_csm_contour_map_ce(wks,vars_3D_diff(ivar,icase,{lev(ilev)},{latS:latE},{lonS:lonE}),res) 
	end do ; do icase=
end do ; do ivar=

plotall	= (/plotdiff(:,0),plotdiff(:,1),plotdiff(:,2)/)
gsn_panel(wks,ndtooned(plotall),(/3,3/),resP)
delete(plotall)


delete(res)
delete(plot)
delete(wks)

end do ; do ilev=

end do ; do is=

 
                               
end
