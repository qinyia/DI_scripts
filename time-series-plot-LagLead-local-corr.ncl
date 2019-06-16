
load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************

;season=(/"ANN","DJF","MAM","JJA","SON"/)
season			=(/"ANN"/)

;compare_type 		= (/"amip-obs","cmip-obs"/)
compare_type 		= (/"amip-obs"/)

compset			= (/"atm"/)

var_TimeLatLon	= (/(/"NetSfcHeat","FSNS","TAUX","LHFLX"/),\
                    (/"TAUX","PSL","NetAtmHeat","LHFLX"/),\
		    (/"CLDLOW","FSNS","TAUX","PSL"/)\
		  /)


do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

	if(compare_type(icomp).eq."amip-obs")then
;		var_TimeLatLon	= (/(/"NetSfcHeat","FLNS","FSNS","LHFLX"/),\
;				    (/"NetAtmHeat","NetSfcHeat","NetToaHeat","LHFLX"/),\
;				    (/"TS","NetAtmHeat","NetSfcHeat","NetToaHeat"/),\
;				    (/"TAUX","PSL","TS","NetSfcHeat"/)\
;				    /)

	end if
	if(compare_type(icomp).eq."cmip-obs")then
;		var_TimeLatLon	= (/(/"NetSfcHeat","FLNS","FSNS","LHFLX"/),\
;				    (/"NetAtmHeat","NetSfcHeat","NetToaHeat","LHFLX"/),\
;				    (/"TS","NetAtmHeat","NetSfcHeat","NetToaHeat"/),\
;				    (/"TAUX","PSL","TS","NetSfcHeat"/)\
;				    /)
	end if

	printVarSummary(var_TimeLatLon)


do icompset=0,dimsizes(compset)-1

;========================================================
; Read data
;========================================================
ncdf = addfile("vars_TimeLatLon_"+season+"_"+compset(icompset)+"_"+compare_type(icomp)+".nc" ,"r")
dSizes				= getfiledimsizes(ncdf)  

do ivarlist=0,dimsizes(var_TimeLatLon(:,0))-1

vars_TimeLatLon		= new((/dimsizes(var_TimeLatLon(ivarlist,:)),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_TimeLatLon!0	= "vars"
vars_TimeLatLon&vars	= ispan(0,dimsizes(var_TimeLatLon(ivarlist,:))-1,1)

do ivar=0,dimsizes(var_TimeLatLon(ivarlist,:))-1
	vars_TimeLatLon(ivar,:,:,:,:)		=ncdf->$var_TimeLatLon(ivarlist,ivar)$ 
end do

printVarSummary(vars_TimeLatLon)
printMinMax(vars_TimeLatLon,True)

;========================================================
; Get regional average of variables
;========================================================

latS				= -20.
latE				= -5.
lonS				= 240.
lonE				= 260.

;latS				= -20.
;latE				= -5.
;lonS				= 180.
;lonE				= 220.

;latS				= -10.
;latE				= -10.
;lonS				= 260.
;lonE				= 260.



lat				= vars_TimeLatLon&lat
lat!0				= "lat"
lat&lat				= vars_TimeLatLon&lat
rad    				= 4.0*atan(1.0)/180.0
clat   				= cos(lat*rad)
copy_VarCoords(lat,clat)

lon				= vars_TimeLatLon&lon

ind_latS			= ind_nearest_coord(latS,lat,0)
ind_latE			= ind_nearest_coord(latE,lat,0)

ind_lonS			= ind_nearest_coord(lonS,lon,0)
ind_lonE			= ind_nearest_coord(lonE,lon,0)

print(ind_latE-ind_latS)
print(ind_lonE-ind_lonS)

icount				= 0


do ilat=ind_latS,ind_latE
do ilon=ind_lonS,ind_lonE


;vars_TimeLatLon_areaave		= wgt_areaave_Wrap(vars_TimeLatLon(:,:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1)

vars_TimeLatLon_areaave		= wgt_areaave_Wrap(vars_TimeLatLon(:,:,:,ilat:ilat,ilon:ilon),clat(ilat:ilat),1.0,1)
printVarSummary(vars_TimeLatLon_areaave)


vars_TimeLatLon_areaave_diff_tmp	= vars_TimeLatLon_areaave(:,1,:)-vars_TimeLatLon_areaave(:,0,:)
;do runave calculation
;vars_TimeLatLon_areaave_diff	= runave_n(vars_TimeLatLon_areaave_diff_tmp,12,0,1)
vars_TimeLatLon_areaave_diff	= vars_TimeLatLon_areaave_diff_tmp

copy_VarCoords_1(vars_TimeLatLon_areaave(vars|:,time|:,case|:),vars_TimeLatLon_areaave_diff)
printVarSummary(vars_TimeLatLon_areaave_diff)

;========================================================
; Do Lead-lag calculations
;========================================================
mxlag				= 12

Index				= dimsizes(vars_TimeLatLon_areaave_diff)
ccr				= new((/Index(0),Index(0),2*mxlag+1/),float)
x_Lead_y			= new((/Index(0),Index(0),mxlag+1/),float)
y_Lead_x			= new((/Index(0),Index(0),mxlag+1/),float)
ccr!0				= "raw_var"
ccr&raw_var			= ispan(0,Index(0)-1,1)
ccr!1				= "other_var"
ccr&other_var			= ispan(0,Index(0)-1,1)
ccr!2				= "time"
ccr&time			= ispan(0,(2*mxlag+1)-1,1)

copy_VarCoords_1(ccr,x_Lead_y)
copy_VarCoords_1(ccr,y_Lead_x)

x_Lead_y!2			= "time"
x_Lead_y&time			= ispan(0,mxlag+1-1,1)
y_Lead_x!2			= "time"
y_Lead_x&time			= ispan(0,mxlag+1-1,1)


do ivar=0,dimsizes(var_TimeLatLon(ivarlist,:))-1
	x_Lead_y(ivar,:,:)			= esccr(vars_TimeLatLon_areaave_diff(ivar,:),vars_TimeLatLon_areaave_diff(:,:),mxlag)
	y_Lead_x(ivar,:,:)			= esccr(vars_TimeLatLon_areaave_diff(:,:),vars_TimeLatLon_areaave_diff(ivar,:),mxlag)

end do

ccr(:,:,0:mxlag-1)			= y_Lead_x(:,:,1:mxlag:-1) ; "negative lag", -1, reverses order
ccr(:,:,mxlag:)				= x_Lead_y(:,:,0:mxlag) ; "positive lag"


;ccr&time			= ispan(0,(2*mxlag+1)-1,1)

printVarSummary(ccr)

; put all data in (latS:latE, lonS:lonE) into one variable: ccr_final, in order to plot in the same figures.
if(ilat.eq.ind_latS.and.ilon.eq.ind_lonS)then
	newdims				= (ind_latE-ind_latS+1)*(ind_lonE-ind_lonS+1)
	ccr_final			= new((/Index(0),newdims,(2*mxlag+1)/),float)
	ccr_final!0			= "raw_var"
	ccr_final&raw_var		= ispan(0,Index(0)-1,1)
	ccr_final!1			= "samples"
	ccr_final&samples		= ispan(0,newdims-1,1)
	ccr_final!2			= "time"
	ccr_final&time			= ispan(0,(2*mxlag+1)-1,1)

end if
printVarSummary(ccr_final)

ccr_final(:,icount,:)			= (/ccr(:,0,:)/) ; (raw_var, other_var, time)

icount					= icount+1

end do ; do ilon=
end do ; do ilat=

printVarSummary(ccr_final)

;========================================================
; Plotting starts here....
;========================================================


wks				= gsn_open_wks("pdf","./figure/time-series-plot_LagLead_Diff_"+compare_type(icomp)+"_"+var_TimeLatLon(ivarlist,0)+"_"+latS+"_"+latE+"_"+lonS+"_"+lonE)

;wks				= gsn_open_wks("x11","./figure/time-series-plot_LagLead_Diff_"+compare_type(icomp)+"_"+var_TimeLatLon(ivarlist,0)+"_"+sprintf("%4.1f",lat(ilat))+"_"+sprintf("%4.1f",lon(ilon)))

plot				= new(4,graphic)
plotd				= new(4,graphic)

;************************************************
; create plots  
;************************************************
trxmin = 0         ; min/max values for X axis 
trxmax = (2*mxlag+1)-1 


xindexc					= ispan(trxmin,trxmax,4)
xlabelc					= xindexc-trxmax/2.

colors				= (/"red","blue","forestgreen","purple"/)

res				= True          ; plot mods desired
res@gsnDraw			= False
res@gsnFrame			= False

res@vpWidthF			= 0.70       ; change aspect ratio of plot
res@vpHeightF			= 0.25

res@tiMainString		="all are correlated to the first var."


nps				= dimsizes(ccr_final&samples)
print(nps)

;res@xyLineColors		= colors
nsteps				= 256/nps
res@xyLineColors	 	=ispan(2,256,nsteps)
res@xyDashPattern		= 0
;res@xyLineThicknessF		= 3.0

res@tiYAxisString		= var_TimeLatLon(ivarlist,0)
res@tiYAxisFontColor		= colors(0)
res@tiYAxisFontHeightF		= 0.05

res@tfPolyDrawOrder		= "PreDraw"    ; Draw any attached primitives
                                         ; before rest of plot.

res@tmXBOn                	= False
res@tmYLLabelFontHeightF  	= 0.02
res@tmYLLabelFontColor    	= res@tiYAxisFontColor

;res@gsnYRefLine			= (/-0.5,0,0.5/)
;res@gsnYRefLineThicknesses	= (/3.0,3.0,3.0/)
;res@gsnYRefLineDashPatterns	= (/2,2,2/)
;
;res@gsnXRefLine			= (/mxlag/)
;res@gsnXRefLineThicknesses	= (/3.0/)
;res@gsnXRefLineDashPatterns	= (/2/)


res@trYMinF			= -1.0
res@trYMaxF			= 1.0

res1				= True

res1@gsnDraw			= False
res1@gsnFrame			= False

res1@vpWidthF			= 0.70
res1@vpHeightF			= 0.25

res1@trYMinF			= -1.0
res1@trYMaxF			= 1.0


res1@xyDashPatterns		= 0                   ; dashed line for 2nd
res1@xyLineColors		= res@xyLineColors
;res1@xyLineColor		= res@xyLineColor
;res1@xyLineThicknessF		= 3.0

res1@tiYAxisFontHeightF		= 0.05
res1@tiYAxisSide		= "Right"	
res1@tiYAxisString		= var_TimeLatLon(ivarlist,1)
res1@tiYAxisFontColor		= colors(1)

res1@tfPolyDrawOrder		= "PreDraw"    ; Draw any attached primitives
                                          ; before rest of plot.

res1@tmXBOn			= False
res1@tmXTOn			= False
res1@tmYLLabelFontColor		= res1@tiYAxisFontColor
res1@tmYLLabelsOn		= False
res1@tmYRLabelFontHeightF	= 0.01      ; resize tick labels
res1@tmYRLabelsOn		= True

;res1@gsnYRefLine		= (/-0.5,0.0,0.5/)
;res1@gsnYRefLineThicknesses	= (/3.0,3.0,3.0/)
;res1@gsnYRefLineDashPatterns	= (/2,2,2/)
;
;res1@gsnXRefLine		= (/mxlag/)
;res1@gsnXRefLineThicknesses	= (/3.0/)
;res1@gsnXRefLineDashPatterns	= (/2/)


;plotd(0)			= gsn_csm_y(wks,ccr(0,1,:),False)

plotd(0)			= gsn_csm_y(wks,ccr_final(0,:,:),res)

;plotd(1)			= gsn_csm_y(wks,ccr(1,:,:),res1)
plotd(1)			= gsn_csm_y(wks,ccr_final(1,:,:),res1)


res@tmXTOn			= False
res@tiYAxisString		= var_TimeLatLon(ivarlist,2)
res@tiYAxisFontColor		= colors(2)
res@tmYLLabelFontColor		= res@tiYAxisFontColor

;plotd(2)			= gsn_csm_y(wks,ccr(2,:,:),res)
plotd(2)			= gsn_csm_y(wks,ccr_final(2,:,:),res)

res1@tmXBOn			= True
res1@tmXBMode			= "Explicit"
res1@tmXBValues			= xindexc  
res1@tmXBLabels			= xlabelc  
res1@tmXBLabelFontHeightF	= 0.02
res1@tiYAxisString		= var_TimeLatLon(ivarlist,3)
res1@tiXAxisString		= "time(month)"
res1@tiYAxisFontColor		= colors(3)
res1@tmYLLabelFontColor		= res1@tiYAxisFontColor

;plotd(3)			= gsn_csm_y(wks,ccr(3,:,:),res1)
plotd(3)			= gsn_csm_y(wks,ccr_final(3,:,:),res1)


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


;end do; do ivar=

end do; do ivarlist=
end do ; do icompset=
end do; do icomp=
                               
end
