
load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************

;season=(/"ANN","DJF","MAM","JJA","SON"/)
season			=(/"ANN"/)

compare_type 		= (/"amip-obs","cmip-obs"/)
;compare_type 		= (/"amip-obs"/)

compset			= (/"atm"/)

var_TimeLatLon	= (/"PRECT","TAUX","PSL","NetAtmHeat","ACRF","LHFLX","CLDLOW","FSNS"/)
;var_TimeLatLon	= (/"CLDLOW","PSL","PRECT","NetAtmHeat","ACRF","LHFLX","TAUX","FSNS"/)


do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

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

latS				= -20.
latE				= -5.
lonS				= 180.
lonE				= 220.

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


vars_TimeLatLonDiff			= vars_TimeLatLon(:,1,:,:,:)-vars_TimeLatLon(:,0,:,:,:)
copy_VarCoords_1(vars_TimeLatLon(vars|:,time|:,lat|:,lon|:,case|:),vars_TimeLatLonDiff)

yvar					= vars_TimeLatLonDiff(vars|:,lat|:,lon|:,time|:)
printVarSummary(yvar)

;do ilat=ind_latS,ind_latE
;do ilon=ind_lonS,ind_lonE


vars_TimeLatLon_areaave		= wgt_areaave_Wrap(vars_TimeLatLon(:,:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1)

;vars_TimeLatLon_areaave		= wgt_areaave_Wrap(vars_TimeLatLon(:,:,:,ilat:ilat,ilon:ilon),clat(ilat:ilat),1.0,1)
printVarSummary(vars_TimeLatLon_areaave)


vars_TimeLatLon_areaave_diff_tmp	= vars_TimeLatLon_areaave(:,1,:)-vars_TimeLatLon_areaave(:,0,:)
;do runave calculation
;vars_TimeLatLon_areaave_diff	= runave_n(vars_TimeLatLon_areaave_diff_tmp,12,0,1)
vars_TimeLatLon_areaave_diff	= vars_TimeLatLon_areaave_diff_tmp

copy_VarCoords_1(vars_TimeLatLon_areaave(vars|:,time|:,case|:),vars_TimeLatLon_areaave_diff)
printVarSummary(vars_TimeLatLon_areaave_diff)

xvar				= vars_TimeLatLon_areaave_diff
printVarSummary(xvar)

;========================================================
; Do Lead-lag calculations
;========================================================
mxlag				= 12

Index				= dimsizes(yvar)
ccr				= new((/Index(0),Index(0),Index(1),Index(2),2*mxlag+1/),float)
x_Lead_y			= new((/Index(0),Index(0),Index(1),Index(2),mxlag+1/),float)
y_Lead_x			= new((/Index(0),Index(0),Index(1),Index(2),mxlag+1/),float)
ccr!0				= "raw_var"
ccr&raw_var			= ispan(0,Index(0)-1,1)
ccr!1				= "other_var"
ccr&other_var			= ispan(0,Index(0)-1,1)
ccr!2				= "lat"
ccr&lat				= yvar&lat
ccr!3				= "lon"
ccr&lon				= yvar&lon
ccr!4				= "time"
ccr&time			= ispan(0,(2*mxlag+1)-1,1)


copy_VarCoords_1(ccr,x_Lead_y)
copy_VarCoords_1(ccr,y_Lead_x)

printVarSummary(ccr)


x_Lead_y!4			= "time"
x_Lead_y&time			= ispan(0,mxlag+1-1,1)
y_Lead_x!4			= "time"
y_Lead_x&time			= ispan(0,mxlag+1-1,1)


do ivar=0,dimsizes(var_TimeLatLon)-1
	x_Lead_y(ivar,:,:,:,:)			= esccr(xvar(ivar,:),yvar(:,:,:,:),mxlag)
	y_Lead_x(ivar,:,:,:,:)			= esccr(xvar(:,:),yvar(ivar,:,:,:),mxlag)

end do

ccr(:,:,:,:,0:mxlag-1)				= y_Lead_x(:,:,:,:,1:mxlag:-1) ; "negative lag", -1, reverses order
ccr(:,:,:,:,mxlag:)				= x_Lead_y(:,:,:,:,0:mxlag) ; "positive lag"

printVarSummary(ccr)
printMinMax(ccr,True)

siglvl			= 0.05
pr			= rtest(x_Lead_y,Index(3),0)
prob			= where(pr.lt.siglvl,pr,yvar@_FillValue)
copy_VarCoords(x_Lead_y,prob)
prob@_FillValue		= yvar@_FillValue

printVarSummary(prob)
printMinMax(prob,True)

;system("/bin/rm -f ccr.nc")
;ncdf = addfile("ccr.nc" ,"c")
;
;ncdf->x_Lead_y		= x_Lead_y



; put all data in (latS:latE, lonS:lonE) into one variable: ccr_final, in order to plot in the same figures.
;if(ilat.eq.ind_latS.and.ilon.eq.ind_lonS)then
;	newdims				= (ind_latE-ind_latS+1)*(ind_lonE-ind_lonS+1)
;	ccr_final			= new((/Index(0),newdims,(2*mxlag+1)/),float)
;	ccr_final!0			= "raw_var"
;	ccr_final&raw_var		= ispan(0,Index(0)-1,1)
;	ccr_final!1			= "samples"
;	ccr_final&samples		= ispan(0,newdims-1,1)
;	ccr_final!2			= "time"
;	ccr_final&time			= ispan(0,(2*mxlag+1)-1,1)
;
;end if
;printVarSummary(ccr_final)
;
;ccr_final(:,icount,:)			= (/ccr(:,1,:)/) ; (raw_var, other_var, time)
;
;icount					= icount+1

;end do ; do ilon=
;end do ; do ilat=

;========================================================
; Plotting starts here....
;========================================================
lat_min 			= -90.
lat_max 			= 90.
lon_min 			= 0.
lon_max 			= 350.



wks				= gsn_open_wks("x11","./figure/time-series-plot_remote_cor_"+compare_type(icomp)+"_"+var_TimeLatLon(0)+"_"+latS+"_"+latE+"_"+lonS+"_"+lonE)


plot				= new(dimsizes(var_TimeLatLon),graphic)

res 				= True
res@gsnFrame 			= False
res@gsnDraw 			= False
res@gsnSpreadColors 		= True

res@cnFillOn 			= True
res@cnLinesOn 			= False
res@cnLineLabelsOn 		= False

;----------------------------------------------------------
; Set Map region
;----------------------------------------------------------
res@mpLandFillColor 		= "gray"
res@mpMinLonF 			= lon_min
res@mpMaxLonF 			= lon_max
res@mpMinLatF 			= lat_min
res@mpMaxLatF 			= lat_max
res@mpCenterLonF 		= 180
;----------------------------------------------------------
; Set Labelbar
;----------------------------------------------------------
res@lbOrientation 		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= -0.01
res@lbLabelStride    		= 2                        ; every other label
;----------------------------------------------------------
; Set X and Y axis; tickmarks
;----------------------------------------------------------
res@tmXBTickSpacingF 		= 40.
res@tmXBLabelFontHeightF  	= 0.014
res@tmYLLabelFontHeightF  	= 0.014

res@cnLevelSelectionMode 	= "ManualLevels" ; manual levels
res@cnMinLevelValF       	= -1.
res@cnMaxLevelValF       	= 1.
res@cnLevelSpacingF      	= 0.1



do ivar=0,dimsizes(var_TimeLatLon)-1
res@gsnCenterString		= var_TimeLatLon(ivar)
plot(ivar)	= gsn_csm_contour_map_ce(wks, x_Lead_y(0,ivar,:,:,0),res) 

;;;; stippling

res3 = True                            ; res2 probability plots
res3@gsnDraw             = False       ; Do not draw plot
res3@gsnFrame            = False       ; Do not advance frome
res3@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
res3@cnMinLevelValF      = 0.00        ; set min contour level
res3@cnMaxLevelValF      = 1.05        ; set max contour level
res3@cnLevelSpacingF     = 0.02        ; set contour spacing
res3@cnInfoLabelOn       = False
res3@cnLinesOn           = False       ; do not draw contour lines
res3@cnLineLabelsOn      = False       ; do not draw contour labels
res3@cnFillScaleF        = 0.6         ; add extra density
                                       ; add cyclic point
res3@gsnLeftString = ""      
res3@gsnRightString = ""

plot2   = gsn_csm_contour(wks,gsn_add_cyclic_point(prob(0,ivar,:,:,0)), res3) 
plot2   = ShadeLtContour(plot2, 0.05, 17) ; shade all areas < 0.07 contour
overlay (plot(ivar), plot2)

end do ; do ivar=

gsn_panel(wks,plot,(/4,2/),False)
exit

end do ; do icompset=
end do; do icomp=
                               
end
