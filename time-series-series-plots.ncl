
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
;var_TimeLatLon			= (/"TS","PSL","CLDLOW","NetSfcHeat","FSNS","TAUX","TAUY"/)
;var_TimeLatLon			= (/"TS","PSL","TAUX","TAUY"/)
;var_TimeLatLon			= (/"TS","PSL","NetSfcHeat","FSNS"/)
var_TimeLatLon			= (/"TS","NetSfcHeat","PSL","TAUX"/)

;compare_type 		= (/"amip-obs"/)
compare_type 		= (/"amip-obs","cmip-obs"/)

compset			= (/"atm"/)


do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))
	if(compare_type(icomp).eq."amip-obs")then
		var_TimeLatLon	= (/"NetSfcHeat","FSNS","PSL","TAUX"/)
	end if
	if(compare_type(icomp).eq."cmip-obs")then
		var_TimeLatLon	= (/"TS","NetSfcHeat","PSL","TAUX"/)
	end if

do icompset=0,dimsizes(compset)-1

;========================================================
; Read data
;========================================================
ncdf = addfile("vars_TimeLatLon_"+season+"_"+compset(icompset)+"_"+compare_type(icomp)+".nc" ,"r")
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

latS				= -10.
latE				= -10.
lonS				= 240.
lonE				= 240.

lat				= vars_TimeLatLon&lat
lat!0				= "lat"
lat&lat				= vars_TimeLatLon&lat
rad    				= 4.0*atan(1.0)/180.0
clat   				= cos(lat*rad)
copy_VarCoords(lat,clat)


vars_TimeLatLon_areaave		= wgt_areaave_Wrap(vars_TimeLatLon(:,:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1)
printVarSummary(vars_TimeLatLon_areaave)


vars_TimeLatLon_areaave_diff_tmp	= vars_TimeLatLon_areaave(:,1,:)-vars_TimeLatLon_areaave(:,0,:)
vars_TimeLatLon_areaave_diff	= runave_n(vars_TimeLatLon_areaave_diff_tmp,6,0,1)

copy_VarCoords_1(vars_TimeLatLon_areaave(vars|:,time|:,case|:),vars_TimeLatLon_areaave_diff)
printVarSummary(vars_TimeLatLon_areaave_diff)

;========================================================
; Plotting starts here....
;========================================================

;do ivar=0,dimsizes(var_TimeLatLon)-1

wks				= gsn_open_wks("pdf","./figure/time-series-series-plots_"+compare_type(icomp)+"_"+latS+"_"+latE+"_"+lonS+"_"+lonE)

plot				= new(4,graphic)
plotd				= new(4,graphic)

;************************************************
; create plots  
;************************************************
trxmin = 0         ; min/max values for X axis 
trxmax = 144

xindexc				= ispan(trxmin,trxmax,12)
xlabelc				= ispan(trxmin,trxmax,12)

res				= True          ; plot mods desired
res@gsnDraw			= False
res@gsnFrame			= False

res@vpWidthF			= 0.70       ; change aspect ratio of plot
res@vpHeightF			= 0.25

res@trXMinF			= trxmin
res@trXMaxF			= trxmax

res@xyLineColor			= "red"
res@xyLineThicknessF		= 2.0

res@tiYAxisString		= var_TimeLatLon(0)
res@tiYAxisFontColor		= res@xyLineColor
res@tiYAxisFontHeightF		= 0.05

res@tfPolyDrawOrder		= "PreDraw"    ; Draw any attached primitives
                                         ; before rest of plot.

res@tmXBOn                	= False
res@tmYLLabelFontHeightF  	= 0.02      ; resize tick labels
res@tmYLLabelFontColor    	= res@xyLineColor

res@gsnYRefLine			= (/0.0/)
res@gsnYRefLineThicknesses	= (/3.0/)
res@gsnYRefLineDashPatterns	= (/2/)



res1				= True

res1@gsnDraw			= False
res1@gsnFrame			= False

res1@vpWidthF			= 0.70
res1@vpHeightF			= 0.25

res1@trXMinF			= trxmin
res1@trXMaxF			= trxmax

res1@xyDashPatterns		= 0                   ; dashed line for 2nd
res1@xyLineColor		= "blue"  ; 3 different colors
res1@xyLineThicknessF		= 2.0

res1@tiYAxisFontHeightF		= 0.05
res1@tiYAxisSide		= "Right"	
res1@tiYAxisString		= var_TimeLatLon(1)
res1@tiYAxisFontColor		= res1@xyLineColor

res1@tfPolyDrawOrder		= "PreDraw"    ; Draw any attached primitives
                                          ; before rest of plot.

res1@tmXBOn			= False
res1@tmXTOn			= False
res1@tmYLLabelFontColor		= res1@xyLineColor
res1@tmYLLabelsOn		= False
res1@tmYRLabelFontHeightF	= 0.01      ; resize tick labels
res1@tmYRLabelsOn		= True

res1@gsnYRefLine		= (/0.0/)
res1@gsnYRefLineThicknesses	= (/3.0/)
res1@gsnYRefLineDashPatterns	= (/2/)

;plot(0) 			= gsn_csm_y(wks,vars_TimeLatLon_areaave(0,:,:),res)
;plotd(0)			= gsn_csm_y(wks,runave_n(vars_TimeLatLon_areaave_diff(0,:),3,0,0),res)
plotd(0)			= gsn_csm_y(wks,vars_TimeLatLon_areaave_diff(0,:),res)

;plot(1) 			= gsn_csm_y(wks,vars_TimeLatLon_areaave(1,:,:),res1)
plotd(1)			= gsn_csm_y(wks,vars_TimeLatLon_areaave_diff(1,:),res1)

res@tmXTOn			= False
res@tiYAxisString		= var_TimeLatLon(2)
res@xyLineColor			= "forestgreen"
res@tiYAxisFontColor		= res@xyLineColor
res@tmYLLabelFontColor		= res@xyLineColor

;plot(2) 			= gsn_csm_y(wks,vars_TimeLatLon_areaave(2,:,:),res)
plotd(2)			= gsn_csm_y(wks,vars_TimeLatLon_areaave_diff(2,:),res)

res1@tmXBOn			= True
res1@tmXBMode			= "Explicit"
res1@tmXBValues			= xindexc  
res1@tmXBLabels			= xlabelc  
res1@tmXBLabelFontHeightF	= 0.02
res1@xyLineColor		= "darkorange3"
res1@tiYAxisString		= var_TimeLatLon(3)
res1@tiXAxisString		= "time(month)"
res1@tiYAxisFontColor		= res1@xyLineColor
res1@tmYLLabelFontColor		= res1@xyLineColor

;plot(3) 			= gsn_csm_y(wks,vars_TimeLatLon_areaave(3,:,:),res1)
plotd(3)			= gsn_csm_y(wks,vars_TimeLatLon_areaave_diff(3,:),res1)

;;;;;; Attach plots along the X axes
attachres1			= True
attachres1@gsnAttachPlotsXAxis	= True  ;; attaches along x-axis
attachres1@gsnAttachBorderOn	= False ;; No border please
attachres2			= True
attachres2@gsnAttachPlotsXAxis	= True  ;; attaches along x-axis
attachres2@gsnAttachBorderOn	= False ;; No border please

;attachid1  			= gsn_attach_plots(plot(0),plot(1:3),attachres1,attachres2)

attachid2  			= gsn_attach_plots(plotd(0),plotd(1:3),attachres1,attachres2)

;;;;;; Maximize output on the page and draw everything
  pres = False               ; No resources needed
  maximize_output(wks,pres)  ; Maximize plot on page


end do ; do icompset=
end do; do icomp=
                               
end
