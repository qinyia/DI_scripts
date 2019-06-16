
load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

;----------------------------------------------------------------------
; Function to attach a labelbar outside of an XY plot
;----------------------------------------------------------------------
function attach_labelbar(wks,plot,labels,colors)
local lbres, vph, vpw, nboxes
begin
  nboxes = dimsizes(labels)
  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False          ; Necessary to control sizes
  lbres@lbPerimOn          = False
  lbres@vpWidthF           = 0.2 * vpw     ; labelbar width
  lbres@vpHeightF          = vph           ; labelbar height
  lbres@lbFillColors       = colors
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.02          ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
  lbres@lbBoxLinesOn       = False

  lbid = gsn_create_labelbar(wks,nboxes,labels,lbres)
;
; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot.
;                                                                 
  amres                = True
  amres@amParallelPosF = 0.61     ; Move away from plot
  annoid = gsn_add_annotation(plot,lbid,amres)

  return(annoid)
  delete([/lbres,vph,vpw,nboxes,lbid,amres,annoid/])
end


;============
; Main codes
;============

begin

season			=(/"ANN"/)

;compare_type 		= (/"amip-obs","cmip-obs"/)
compare_type 		= (/"amip-obs"/)

compset			= (/"atm"/)

;var_TimeLatLon	= (/"NetSfcHeat","TAUX","PSL","LHFLX"/)
;var_TimeLatLon	= (/"NetSfcHeat","LHFLX"/)
var_TimeLatLon	= (/"NetAtmHeat","PRECT"/)


do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))


do icompset=0,dimsizes(compset)-1

;========================================================
; Read data
;========================================================
ncdf = addfile("vars_TimeLatLon_"+season+"_"+compset(icompset)+"_"+compare_type(icomp)+".nc" ,"r")
dSizes				= getfiledimsizes(ncdf)  

do ivarlist=0,dimsizes(var_TimeLatLon)-1

vars_TimeLatLon		= new((/dimsizes(var_TimeLatLon),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_TimeLatLon!0	= "vars"
vars_TimeLatLon&vars	= ispan(0,dimsizes(var_TimeLatLon)-1,1)

do ivar=0,dimsizes(var_TimeLatLon)-1
	vars_TimeLatLon(ivar,:,:,:,:)		=ncdf->$var_TimeLatLon(ivar)$ 
end do

end do; do ivarlist=

printVarSummary(vars_TimeLatLon)
printMinMax(vars_TimeLatLon,True)

;========================================================
; Get regional average of variables
;========================================================

latS				= -20.
latE				= -10.
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

time				= vars_TimeLatLon&time

ind_latS			= ind_nearest_coord(latS,lat,0)
ind_latE			= ind_nearest_coord(latE,lat,0)

ind_lonS			= ind_nearest_coord(lonS,lon,0)
ind_lonE			= ind_nearest_coord(lonE,lon,0)

print(ind_latE-ind_latS)
print(ind_lonE-ind_lonS)

icount				= 0


;do ilat=ind_latS,ind_latE
;do ilon=ind_lonS,ind_lonE


;---------------------------------------------------------------
; get variable from the specific region
;---------------------------------------------------------------
; vars_TimeLatLon(var,case,time,lat,lon)
vars_TimeLatLon_areaave					= vars_TimeLatLon(:,:,{72::3},ind_latS:ind_latE,ind_lonS:ind_lonE)
printVarSummary(vars_TimeLatLon_areaave)
;---------------------------------------------------------------
; get the difference
;---------------------------------------------------------------
vars_TimeLatLon_areaave_diff				= vars_TimeLatLon_areaave(:,1,:,:,:)-vars_TimeLatLon_areaave(:,0,:,:,:)
copy_VarCoords_1(vars_TimeLatLon_areaave(vars|:,time|:,lat|:,lon|:,case|:),vars_TimeLatLon_areaave_diff)
printVarSummary(vars_TimeLatLon_areaave_diff)
;---------------------------------------------------------------
; get index
;---------------------------------------------------------------
Index							= dimsizes(vars_TimeLatLon_areaave)
printVarSummary(Index)
print(Index)
;---------------------------------------------------------------
; (vars,case,time*lat*lon)
;---------------------------------------------------------------
newdims1						= dimsizes(time(72::3))*(ind_latE-ind_latS+1)*(ind_lonE-ind_lonS+1)
vars_TimeLatLon_areaave_LatLon				= reshape(vars_TimeLatLon_areaave,(/Index(0),Index(1),newdims1/))
vars_TimeLatLon_areaave_LatLon!0			= "vars"
vars_TimeLatLon_areaave_LatLon&vars			= ispan(0,Index(0)-1,1)
vars_TimeLatLon_areaave_LatLon!1			= "case"
vars_TimeLatLon_areaave_LatLon&case			= ispan(0,Index(1)-1,1)
vars_TimeLatLon_areaave_LatLon!2			= "TimeLatLon"
vars_TimeLatLon_areaave_LatLon&TimeLatLon		= ispan(0,newdims1-1,1)
printVarSummary(vars_TimeLatLon_areaave_LatLon)

vars_TimeLatLon_areaave_diff_LatLon			= reshape(vars_TimeLatLon_areaave_diff,(/Index(0),newdims1/))
copy_VarCoords_1(vars_TimeLatLon_areaave_LatLon(vars|:,TimeLatLon|:,case|:),vars_TimeLatLon_areaave_diff_LatLon)
printVarSummary(vars_TimeLatLon_areaave_diff_LatLon)
;---------------------------------------------------------------
; (vars,case,time*lat,lon)
;---------------------------------------------------------------
newdims2						= dimsizes(time(72::3))*(ind_latE-ind_latS+1)
vars_TimeLatLon_areaave_Lat				= reshape(vars_TimeLatLon_areaave,(/Index(0),Index(1),newdims2,Index(4)/))

vars_TimeLatLon_areaave_Lat!0				= "vars"
vars_TimeLatLon_areaave_Lat&vars			= ispan(0,Index(0)-1,1)
vars_TimeLatLon_areaave_Lat!1				= "case"
vars_TimeLatLon_areaave_Lat&case			= ispan(0,Index(1)-1,1)
vars_TimeLatLon_areaave_Lat!2				= "TimeLat"
vars_TimeLatLon_areaave_Lat&TimeLat			= ispan(0,newdims2-1,1)
vars_TimeLatLon_areaave_Lat!3				= "lon"
vars_TimeLatLon_areaave_Lat&lon				= lon(ind_lonS:ind_lonE)
printVarSummary(vars_TimeLatLon_areaave_Lat)

vars_TimeLatLon_areaave_diff_Lat			= reshape(vars_TimeLatLon_areaave_diff,(/Index(0),newdims2,Index(4)/))
copy_VarCoords_1(vars_TimeLatLon_areaave_Lat(vars|:,TimeLat|:,lon|:,case|:),vars_TimeLatLon_areaave_diff_Lat)
printVarSummary(vars_TimeLatLon_areaave_diff_Lat)
;---------------------------------------------------------------
; (vars,case,time*lon,lat)
;---------------------------------------------------------------
newdims3						= dimsizes(time(72::3))*(ind_lonE-ind_lonS+1)
vars_TimeLatLon_areaave_Lon				= reshape(vars_TimeLatLon_areaave,(/Index(0),Index(1),newdims3,Index(3)/))

vars_TimeLatLon_areaave_Lon!0				= "vars"
vars_TimeLatLon_areaave_Lon&vars			= ispan(0,Index(0)-1,1)
vars_TimeLatLon_areaave_Lon!1				= "case"
vars_TimeLatLon_areaave_Lon&case			= ispan(0,Index(1)-1,1)
vars_TimeLatLon_areaave_Lon!2				= "TimeLon"
vars_TimeLatLon_areaave_Lon&TimeLon			= ispan(0,newdims3-1,1)
vars_TimeLatLon_areaave_Lon!3				= "lat"
vars_TimeLatLon_areaave_Lon&lat				= lat(ind_latS:ind_latE)
printVarSummary(vars_TimeLatLon_areaave_Lon)

vars_TimeLatLon_areaave_diff_Lon			= reshape(vars_TimeLatLon_areaave_diff,(/Index(0),newdims3,Index(3)/))
copy_VarCoords_1(vars_TimeLatLon_areaave_Lon(vars|:,TimeLon|:,lat|:,case|:),vars_TimeLatLon_areaave_diff_Lon)
printVarSummary(vars_TimeLatLon_areaave_diff_Lon)

;---------------------------------------------------------------
; Set the plotting variables
;---------------------------------------------------------------
xx							= vars_TimeLatLon_areaave_Lon(0,0,:,:)
yy							= vars_TimeLatLon_areaave_Lon(1,0,:,:)
printVarSummary(xx)

levels							= xx&lat
labels							= sprintf("%3.1f",levels)
print(labels)

inds							= dimsizes(xx)
nps							= inds(1)

nsteps							= 256/nps
colors							= ispan(2,256,nsteps)
overlays 						= new(nps-1,graphic)
overlaysDiff 						= new(nps-1,graphic)

xxDiff							= vars_TimeLatLon_areaave_diff_Lon(0,:,:)
yyDiff							= vars_TimeLatLon_areaave_diff_Lon(1,:,:)
printVarSummary(xxDiff)

;========================================================
; Plotting starts here....
;========================================================
wks				= gsn_open_wks("x11","./figure/scatter_"+compare_type(icomp)+"_"+var_TimeLatLon(0)+"_"+var_TimeLatLon(1)+"_"+latS+"_"+latE+"_"+lonS+"_"+lonE)


;************************************************
; create plots  
;************************************************


res				= True

res@gsnFrame			= False
res@gsnDraw			= False

res@xyMarkLineModes		= "Markers" 
res@xyMarkers			=  16 
res@xyMarkerSizeF		= 0.01 
res@tmLabelAutoStride 		= True 

res@xyMarkerColor		= colors(0)

res@tiXAxisString		= var_TimeLatLon(0)
res@tiYAxisString		= var_TimeLatLon(1)


;==========================================
; Difference plots
;==========================================
res@trYMinF			= min(yyDiff)
res@trYMaxF			= max(yyDiff)
res@trXMinF			= min(xxDiff)
res@trXMaxF			= max(xxDiff)

plotd				= gsn_csm_xy(wks,xxDiff(:,0),yyDiff(:,0),res)

do iloop = 1,nps-1
	res@xyMarkerColor		= colors(iloop)
	overlaysDiff(iloop-1)		= gsn_csm_xy(wks,xxDiff(:,iloop),yyDiff(:,iloop),res)

	overlay(plotd,overlaysDiff(iloop-1))
end do

labelbarDiff			= attach_labelbar(wks,plotd,labels,colors)

maximize_output(wks,False)

;==========================================
; Climatology plots
;==========================================
res@trYMinF			= min(yy)
res@trYMaxF			= max(yy)
res@trXMinF			= min(xx)
res@trXMaxF			= max(xx)

plot				= gsn_csm_xy(wks,xx(:,0),yy(:,0),res)

do iloop = 1,nps-1
	res@xyMarkerColor		= colors(iloop)
	overlays(iloop-1)		= gsn_csm_xy(wks,xx(:,iloop),yy(:,iloop),res)

	overlay(plot,overlays(iloop-1))
end do

labelbar 			= attach_labelbar(wks,plot,labels,colors)

maximize_output(wks,False)



exit


end do ; do icompset=
end do; do icomp=
                               
end
