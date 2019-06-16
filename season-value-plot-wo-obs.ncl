
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

begin


var_SeaLatLon=(/"FSNS","TS","PSL","TAUX"/) 
;compare_type = (/"amip-obs","cmip-obs","cmip-amip"/)
;compare_type = (/"cmip-amip"/)
compare_type = (/"cmip-obs"/)

case=(/"def","new"/)

Dir="/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

ncdf 			= addfile("vars_SeaLatLon_"+compare_type(icomp)+"-wo-obs.nc","r")

dSizes			= getfiledimsizes(ncdf)

vars_SeaLatLon		= new((/dimsizes(var_SeaLatLon),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_SeaLatLon!0	= "vars"
vars_SeaLatLon&vars	= ispan(0,dimsizes(var_SeaLatLon)-1,1)

do ivar=0,dimsizes(var_SeaLatLon)-1
	vars_SeaLatLon(ivar,:,:,:,:)		= ncdf->$var_SeaLatLon(ivar)$
end do

printVarSummary(vars_SeaLatLon)

;========================================================================
; Set plotting parameters
;========================================================================

latS 				= -20
latE 				= -10

lonS				= 220
lonE				= 260

TimeLat				= False
series_plots			= True

;========================================================================
; Area Average mean over specific region
;========================================================================
lat					= vars_SeaLatLon&lat
lat!0					= "lat"
lat&lat					= vars_SeaLatLon&lat

rad    					= 4.0*atan(1.0)/180.0
clat   					= cos(lat*rad)
copy_VarCoords(lat,clat)

vars_SeaLatLon_areaave 			= wgt_areaave_Wrap(vars_SeaLatLon(:,:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1)

vars_SeaLatLon_diff_areaave		= vars_SeaLatLon_areaave(:,0,:)
vars_SeaLatLon_diff_areaave		= vars_SeaLatLon_areaave(:,1,:)-vars_SeaLatLon_areaave(:,0,:) ; def-obs

vars_SeaLatLon_areaave_rev 		= vars_SeaLatLon_areaave(vars|:,case|:,time|:)
vars_SeaLatLon_diff_areaave_rev 	= vars_SeaLatLon_diff_areaave(vars|:,time|:)

printVarSummary(vars_SeaLatLon_diff_areaave_rev)
printVarSummary(vars_SeaLatLon_areaave_rev)

delete([/lat,clat/])

;========================================================================
; Plotting starts here....
;========================================================================
wks   = gsn_open_wks ("x11", "./figure/season-value_"+compare_type(icomp)) 


if(series_plots)then
	plot				=new(4,graphic)
	plotd				=new(4,graphic)

	trxmin 				= 0 
	trxmax 				= 11
	
	xindexc				= ispan(trxmin,trxmax,12)
	xlabelc				= ispan(trxmin,trxmax,12)
	
	res				= True
	res@gsnDraw			= False
	res@gsnFrame			= False
	
	res@vpWidthF			= 0.70
	res@vpHeightF			= 0.25

	res@trXMinF			= trxmin
	res@trXMaxF			= trxmax
	
	res@xyLineColor			= "red"
	res@xyLineThicknessF		= 3.0
	res@xyDashPatterns		= (/0,1,2/)
	
	res@tiYAxisString		= var_SeaLatLon(0)
	res@tiYAxisFontColor		= res@xyLineColor
	res@tiYAxisFontHeightF		= 0.05
	
	res@tfPolyDrawOrder		= "PreDraw"
	
	res@tmXBOn                	= False
	res@tmYLLabelFontHeightF  	= 0.01
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
	
;	res1@xyDashPatterns		= 0
	res1@xyLineColor		= "blue"
	res1@xyLineThicknessF		= 3.0
	res1@xyDashPatterns		= (/0,1,2/)

	res1@tiYAxisFontHeightF		= 0.05
	res1@tiYAxisSide		= "Right"	
	res1@tiYAxisString		= var_SeaLatLon(1)
	res1@tiYAxisFontColor		= res1@xyLineColor
	
	res1@tfPolyDrawOrder		= "PreDraw"	

	res1@tmXBOn			= False
	res1@tmXTOn			= False
	res1@tmYLLabelFontColor		= res1@xyLineColor
	res1@tmYLLabelsOn		= False
	res1@tmYRLabelFontHeightF	= 0.01
	res1@tmYRLabelsOn		= True
	
	res1@gsnYRefLine		= (/0.0/)
	res1@gsnYRefLineThicknesses	= (/3.0/)
	res1@gsnYRefLineDashPatterns	= (/2/)
	
	;plot(0) 			= gsn_csm_y(wks,vars_TimeLatLon_areaave(0,:,:),res)
	plotd(0)			= gsn_csm_y(wks,vars_SeaLatLon_diff_areaave(0,:),res)

	;plot(1) 			= gsn_csm_y(wks,vars_TimeLatLon_areaave(1,:,:),res1)
	plotd(1)			= gsn_csm_y(wks,vars_SeaLatLon_diff_areaave(1,:),res1)

	res@tmXTOn			= False
	res@tiYAxisString		= var_SeaLatLon(2)
	res@xyLineColor			= "forestgreen"
	res@tiYAxisFontColor		= res@xyLineColor
	res@tmYLLabelFontColor		= res@xyLineColor
	
	;plot(2) 			= gsn_csm_y(wks,vars_TimeLatLon_areaave(2,:,:),res)
	plotd(2)			= gsn_csm_y(wks,vars_SeaLatLon_diff_areaave(2,:),res)

	res1@tmXBOn			= True
	res1@tmXBMode			= "Explicit"
	res1@tmXBValues			= xindexc  
	res1@tmXBLabels			= xlabelc  
	res1@tmXBLabelFontHeightF	= 0.02
	res1@xyLineColor		= "darkorange3"
	res1@tiYAxisString		= var_SeaLatLon(3)
	res1@tiXAxisString		= "time(month)"
	res1@tiYAxisFontColor		= res1@xyLineColor
	res1@tmYLLabelFontColor		= res1@xyLineColor
	
	;plot(3) 			= gsn_csm_y(wks,vars_TimeLatLon_areaave(3,:,:),res1)
	plotd(3)			= gsn_csm_y(wks,vars_SeaLatLon_diff_areaave(3,:),res1)

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

end if ; series_plots

end do ; do icomp=
                              
end
