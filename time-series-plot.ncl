
load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************

;season=(/"ANN","DJF","MAM","JJA","SON"/)
season			=(/"ANN"/)
;var_TimeLatLon			=(/"TEMP","UVEL","VVEL","WVEL"/)
;var_TimeLatLon			=(/"TEMP"/)
;var_TimeLatLon			= (/"HMXL"/)

;var_TimeLatLon			= (/"PSL","CLDLOW"/)
;var_TimeLatLon			= (/"TS"/)
var_TimeLatLon			= (/"TS","PSL","CLDLOW","NetSfcHeat","FSNS","TAUX","TAUY"/)

compare_type 		= (/"cmip-obs"/)
case			= (/"DEF","NEW"/)

compset			= (/"atm"/)


do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

do icompset=0,dimsizes(compset)-1

;========================================================
; Read data
;========================================================
ncdf = addfile("vars_TimeLatLon_"+season+"_"+compset(icompset)+".nc" ,"r")
dSizes				= getfiledimsizes(ncdf)  

vars_TimeLatLon		= new((/dimsizes(var_TimeLatLon),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_TimeLatLon!0	= "vars"
vars_TimeLatLon&vars	= ispan(0,dimsizes(var_TimeLatLon)-1,1)

do ivar=0,dimsizes(var_TimeLatLon)-1
	vars_TimeLatLon(ivar,:,:,:,:)		=ncdf->$var_TimeLatLon(ivar)$ 
end do

printVarSummary(vars_TimeLatLon)
printMinMax(vars_TimeLatLon,True)

;========================================================
; Get regional average of variables
;========================================================

latS				= -20.
latE				= -5.
lonS				= 180.
lonE				= 220.

lat				= vars_TimeLatLon&lat
lat!0				= "lat"
lat&lat				= vars_TimeLatLon&lat
rad    				= 4.0*atan(1.0)/180.0
clat   				= cos(lat*rad)
copy_VarCoords(lat,clat)


vars_TimeLatLon_areaave		= wgt_areaave_Wrap(vars_TimeLatLon(:,:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1)
printVarSummary(vars_TimeLatLon_areaave)

vars_TimeLatLon_areaave_diff	= vars_TimeLatLon_areaave(:,1,:)-vars_TimeLatLon_areaave(:,0,:)
copy_VarCoords_1(vars_TimeLatLon_areaave(vars|:,time|:,case|:),vars_TimeLatLon_areaave_diff)
printVarSummary(vars_TimeLatLon_areaave_diff)

;========================================================
; Plotting starts here....
;========================================================

do ivar=0,dimsizes(var_TimeLatLon)-1

wks				= gsn_open_wks("x11","./figure/time-evolution")

plot				= new(4,graphic)

res 				= True

res@gsnFrame		 	= False
res@gsnDraw			= False

;----------------------------------------------------
; Plot Height and Width
;----------------------------------------------------
res@vpWidthF			= 0.7
res@vpHeightF			= 0.4
;----------------------------------------------------
; Left/Center/Right Strings
;----------------------------------------------------
res@gsnLeftString		= var_TimeLatLon(ivar)
if(latS.lt.0)then
	if(latE.lt.0)then
		if(lonS.le.180)then
			if(lonE.le.180)then
				res@gsnRightString		= abs(latS)+"S-"+abs(latE)+"S;"+abs(lonS)+"E-"+abs(lonE)+"E"
			else
				res@gsnRightString		= abs(latS)+"S-"+abs(latE)+"S;"+abs(lonS)+"E-"+abs(360-lonE)+"S"
			end if
		else
				if(lonE.le.180)then
				res@gsnRightString		= abs(latS)+"S-"+abs(latE)+"S;"+abs(360-lonS)+"S-"+abs(lonE)+"E"
			else
				res@gsnRightString		= abs(latS)+"S-"+abs(latE)+"S;"+abs(360-lonS)+"S-"+abs(360-lonE)+"S"
			end if
		end if
	else
		if(lonS.le.180)then
			if(lonE.le.180)then
				res@gsnRightString		= abs(latS)+"S-"+abs(latE)+"N;"+abs(lonS)+"E-"+abs(lonE)+"E"
			else
				res@gsnRightString		= abs(latS)+"S-"+abs(latE)+"N;"+abs(lonS)+"E-"+abs(360-lonE)+"S"
			end if
		else
			if(lonE.le.180)then
				res@gsnRightString		= abs(latS)+"S-"+abs(latE)+"N;"+abs(360-lonS)+"S-"+abs(lonE)+"E"
			else
				res@gsnRightString		= abs(latS)+"S-"+abs(latE)+"N;"+abs(360-lonS)+"S-"+abs(360-lonE)+"S"
			end if
		end if

	end if

else
	if(latE.gt.0)then
		if(lonS.le.180)then
			if(lonE.le.180)then
				res@gsnRightString		= abs(latS)+"N - "+abs(latE)+"N;"+abs(lonS)+"E-"+abs(lonE)+"E"
			else
				res@gsnRightString		= abs(latS)+"N - "+abs(latE)+"N;"+abs(lonS)+"E-"+abs(360-lonE)+"S"
			end if
		else
			if(lonE.le.180)then
				res@gsnRightString		= abs(latS)+"N - "+abs(latE)+"N;"+abs(360-lonS)+"S-"+abs(lonE)+"E"
			else
				res@gsnRightString		= abs(latS)+"N - "+abs(latE)+"N;"+abs(360-lonS)+"S-"+abs(360-lonE)+"S"
			end if
		end if

	end if
end if
;----------------------------------------------------
; Turn on/off grid lines
;----------------------------------------------------
;res@tmXMajorGrid                = True          ; implement x grid 
;res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
;res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
;res@tmYMajorGrid                = True          ; implement y grid
;res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
;res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
;----------------------------------------------------
; X/Y Axises and ranges
;----------------------------------------------------
res@trYMinF			= min(vars_TimeLatLon_areaave(ivar,:,:))
res@trYMaxF			= max(vars_TimeLatLon_areaave(ivar,:,:))

xsize				= dimsizes(vars_TimeLatLon_areaave(ivar,0,:))
res@tmXBMode   			= "Explicit"
res@tmXBValues 			= fspan(0,xsize,13)
res@tmXBLabels 			= fspan(0,xsize,13)

res@tiYAxisString		= var_TimeLatLon(ivar)
res@tiXAxisString		= "time(month)"
;----------------------------------------------------
; Line features
;----------------------------------------------------
res@xyLineThicknessF 		= 3.0

res@xyDashPattern		= 0
res@xyLineColors		= (/"red","blue"/)

res@gsnCenterString		= case(0)
plot(0)			 	= gsn_csm_y(wks,vars_TimeLatLon_areaave(ivar,0,:),res)

res@gsnCenterString		= case(1)
plot(1)			 	= gsn_csm_y(wks,vars_TimeLatLon_areaave(ivar,1,:),res)

res@gsnCenterString		= ""
plot(2)			 	= gsn_csm_y(wks,vars_TimeLatLon_areaave(ivar,:,:),res)

;----------------------------------------------------
; X/Y Axises and ranges
;----------------------------------------------------
res@trYMinF			= min(vars_TimeLatLon_areaave_diff(ivar,:))
res@trYMaxF			= max(vars_TimeLatLon_areaave_diff(ivar,:))

xsize				= dimsizes(vars_TimeLatLon_areaave_diff(ivar,:))
res@tmXBMode   			= "Explicit"
res@tmXBValues 			= fspan(0,xsize,13)
res@tmXBLabels 			= fspan(0,xsize,13)

res@gsnCenterString		= "NEW-DEF"
plot(3)			 	= gsn_csm_y(wks,vars_TimeLatLon_areaave_diff(ivar,:),res)

;========start to add legends
gres 				= True
gres@YPosPercent 		= 95.    ; expressed as %, 0->100, sets position of top border of legend 
                          		 ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
gres@XPosPercent 		= 5      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

lineres 			= True
lineres@lgLineColors 		= res@xyLineColors
lineres@lgLineThicknesses 	= res@xyLineThicknessF
lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line


textres 			= True
textres@lgLabels 		= (/"DEF","NEW"/)

plot(2)				= simple_legend(wks,plot(2),gres,lineres,textres)


gsn_panel(wks,plot,(/2,2/),False)


end do; do ivar=
end do ; do icompset=
end do; do icomp=
                               
end
