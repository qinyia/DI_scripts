
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

; 3D
function regrid_diff_3D(var_def[*][*][*]:numeric,lev_def[*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		        var_new[*][*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

var_interp 				= lev_lat_lon_interpolate_2vars(var_def,lev_def,lat_def,lon_def, \
				       				        var_new,lat_new,lon_new)

var_interp_diff 			= var_interp(0,:,:,:)
var_interp_diff 			= var_interp(1,:,:,:) - var_interp(0,:,:,:) ; def-obs

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end

; 2D
function regrid_diff_2D(var_def[*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		        var_new[*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

;horizontal interpolation
var_interp 				= lat_lon_interpolate_2vars(var_def,lat_def,lon_def, \
				       				    var_new,lat_new,lon_new)

var_interp_diff 			= var_interp(0,:,:)
var_interp_diff 			= var_interp(1,:,:) - var_interp(0,:,:) ; def-obs

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin

var2D				= (/"LHFLX_revised","SHFLX","RadCooling_allsky","PS","OCNFRAC"/)

var3D				= (/"U","V","OMEGA","MSE"/)


;compare_type 			= (/"cmip-obs"/)
;compare_type 			= (/"amip-obs"/)
compare_type 			= (/"amip-obs","cmip-obs"/)


;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Setting...
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DataProcess 			= False
LatLonPlot	 		= False
LonValuePlot			= False
AreaMeanPlot			= True
plot_type			= "x11"
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dir				= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"


if(DataProcess)then

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))
	if(compare_type(icomp) .eq."amip-obs")then
		casename = (/"FC5_f09f09_MG15_COSP","new-cldfrac16-1d-MG15-revise-2-COSP"/)
		fdir = (/Dir+"amip/",Dir+"amip/"/)
	end if
	if(compare_type(icomp) .eq."cmip-obs")then
		casename = (/"BC5_f19g16","BC5_f19g16_mac2"/)
		fdir = (/Dir+"atm/",Dir+"atm/"/)
	end if

	fapend 				= (/"_ANN_climo.nc", \
			                    "_ANN_climo.nc" /)

;=====================================================================
; deal with 2D variables
;=====================================================================

do ivar=0,dimsizes(var2D)-1
	print(var2D(ivar))

	result_def				= get_var_2D(fdir(0),fapend(0),casename(0),var2D(ivar),"model")
	var_def					= result_def[0]
	lat_def					= result_def[1]
	lon_def					= result_def[2]
	
	result_new				= get_var_2D(fdir(1),fapend(1),casename(1),var2D(ivar),"model")
	var_new					= result_new[0]
	lat_new					= result_new[1]
	lon_new					= result_new[2]

	delete([/result_def,result_new/])


	result 					= regrid_diff_2D(var_def,lat_def,lon_def,\
		              		                         var_new,lat_new,lon_new)


	var_2D				= result[0]
	var_2D_diff 			= result[1]

	delete([/result,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])

	; put all 2D variables into one variable

	index				= dimsizes(var_2D)
	if(ivar.eq.0)then
		vars_LatLon		= new((/dimsizes(var2D),index(0),index(1),index(2)/),float)
		vars_LatLon!0		= "vars"
		vars_LatLon&vars	= ispan(0,dimsizes(var2D)-1,1)
	end if
	delete(index)

	vars_LatLon(ivar,:,:,:)		= var_2D

	delete([/var_2D,var_2D_diff/])
	
end do ; do ivar=

printVarSummary(vars_LatLon)
;=====================================================================
; deal with 3D variables
;=====================================================================

do ivar=0,dimsizes(var3D)-1
	print(var3D(ivar))

	result_def				= get_var_3D(fdir(0),fapend(0),casename(0),var3D(ivar),"model")
	var_def					= result_def[0]
	lat_def					= result_def[1]
	lon_def					= result_def[2]
	lev_def					= result_def[3]
	
	result_new				= get_var_3D(fdir(1),fapend(1),casename(1),var3D(ivar),"model")
	var_new					= result_new[0]
	lat_new					= result_new[1]
	lon_new					= result_new[2]

	delete([/result_def,result_new/])

	result 					= regrid_diff_3D(var_def,lev_def,lat_def,lon_def,\
		              		                         var_new,lat_new,lon_new)

	var_3D 				= result[0]
	var_3D_diff 			= result[1]

	delete([/result,var_def,lat_def,lon_def,lev_def,var_new,lat_new,lon_new/])

	; put all 2D variables into one variable

	index				= dimsizes(var_3D)
	if(ivar.eq.0)then
		vars_LevLatLon		= new((/dimsizes(var3D),index(0),index(1),index(2),index(3)/),float)
		vars_LevLatLon!0		= "vars"
		vars_LevLatLon&vars	= ispan(0,dimsizes(var3D)-1,1)
	end if
	delete(index)

	vars_LevLatLon(ivar,:,:,:,:)		= var_3D
	
	delete([/var_3D,var_3D_diff/])
end do; do ivar=

printVarSummary(vars_LevLatLon)

;=====================================================================
; initially test  the basic feature of these variables
;=====================================================================
if(False)then
wks			= gsn_open_wks(plot_type,"./figure/"+var2D(0))

cmap			= read_colormap_file("BlueYellowRed")

;gsn_define_colormap(wks,"BlueYellowRed")
;gsn_define_colormap(wks,"MPL_coolwarm")

res 			= True

res@cnFillPalette	= cmap(20:233,:)
res@gsnDraw		= False
res@gsnFrame		= False

res@cnFillOn		= True

;res@cnLevelSelectionMode	= "ManualLevels"
;res@cnMinLevelValF		= -140.
;res@cnMaxLevelValF		= 140.
;res@cnLevelSpacingF		= 20.

plot2			= new(dimsizes(var2D),graphic)
res@mpCenterLonF	= 180.

do ivar2 = 0,dimsizes(var2D)-1
	res@gsnCenterString			= var2D(ivar2)

	maxval                          = max(vars_LatLon(ivar2,0,:,:))
        minval                          = min(vars_LatLon(ivar2,0,:,:))
        spcval                          = (maxval-minval)/11

	tmp                             = (/abs(maxval),abs(minval)/)
	exmval                          = max(tmp)
	data_levels                     = fspan(-1.*exmval,exmval,11)
	delete(tmp)
    
        res@cnLevelSelectionMode        = "ExplicitLevels"
        res@cnLevels                    = data_levels

        res@cnLineLabelFormat           = "0@;*.3f"                      ; Only one value after decimal point
        if(abs(spcval).lt.1.e-4)then
                res@lbLabelStrings              = sprintf("%5.1g",data_levels)   ; Format the labelbar labels
        else
                res@lbLabelStrings              = sprintf("%5.2f",data_levels)   ; Format the labelbar labels
        end if

	plot2(ivar2)			= gsn_csm_contour_map_ce(wks,vars_LatLon(ivar2,0,:,:),res)
end do
gsn_panel(wks,plot2,(/3,2/),False)

exit

end if
;=====================================================================
; Calculating MSE budget terms starts here...
;=====================================================================

latS				= -30.
latE				= 30.
lonS				= 120.
lonE				= 290.

;---------------------------------------------------------------------
; read all necessary data
;---------------------------------------------------------------------
LHFLX				= vars_LatLon(0,:,:,:)
SHFLX				= vars_LatLon(1,:,:,:)
RadCooling_allsky		= vars_LatLon(2,:,:,:)
PS				= vars_LatLon(3,:,:,:)
OCNFRAC				= vars_LatLon(4,:,:,:)


uwind				= vars_LevLatLon(0,:,:,:,:)
vwind				= vars_LevLatLon(1,:,:,:,:)
wwind				= vars_LevLatLon(2,:,:,:,:)
MSE				= vars_LevLatLon(3,:,:,:,:)


lat				= vars_LevLatLon&lat
lon				= vars_LevLatLon&lon
lev				= vars_LevLatLon&lev

delete([/vars_LatLon,vars_LevLatLon/])
;---------------------------------------------------------------------
; external forcing = latent heat + sensible heat + radiative cooling
;---------------------------------------------------------------------
exforcing			= LHFLX+SHFLX+RadCooling_allsky
copy_VarCoords(LHFLX,exforcing)


; test whether external forcing is global-mean conserved: Yes. It is conserved.
if(False)then
	lat_here 				= lat
	rad 					= 4.0*atan(1.0)/180.
	clat 					= cos(lat_here*rad)
	clat!0 					= "lat"
	clat&lat 				= lat_here
	
	globalmean				=  wgt_areaave_Wrap(exforcing,clat,1.0,1)
	print(globalmean)
end if

;---------------------------------------------------------------------
; horizontal advection
;---------------------------------------------------------------------
gradLatLon			= grad_latlon_cfd(MSE,lat,lon,False,False)
DhDy				= gradLatLon[0]
DhDx				= gradLatLon[1]


u_DhDx				= uwind*DhDx*(-1.)
v_DhDy				= vwind*DhDy*(-1.)

copy_VarCoords(uwind,u_DhDx)
copy_VarCoords(vwind,v_DhDy)
printVarSummary(u_DhDx)

;---------------------------------------------------------------------
; vertical  advection
;---------------------------------------------------------------------
DhDp  				= center_finite_diff_n (MSE,lev*100.,False,0,1) ; (kJ/kg)/(Pa)
copy_VarCoords(wwind,DhDp)

w_DhDp				= wwind*DhDp*(-1.)
copy_VarCoords(wwind,w_DhDp)

printVarSummary(w_DhDp)
printMinMax(w_DhDp,True)

;------------------------------------------------------------------
; use plotting to test the spatial pattern
;------------------------------------------------------------------
if(False)then

wks = gsn_open_wks(plot_type,"./figure/test-DhDp")
res = True
res@cnFillOn = True
res@cnLinesOn	= False
res@cnLineLabelsOn = False
res@gsnFrame = False
res@gsnDraw  = False

res@mpMinLonF	= 120.
res@mpMaxLonF	= 290.
res@mpMinLatF	= -30.
res@mpMaxLatF	= 30.

res@mpCenterLonF	 = 180.
res@gsnAddCyclic	= False

plot = new(30,graphic)

do ilev=0,dimsizes(lev)-1
	res@gsnCenterString	= lev(ilev)
	plot(ilev) = gsn_csm_contour_map_ce(wks,DhDp(0,ilev,{-30:30},{120:290}),res)
end do
gsn_panel(wks,plot,(/5,6/),False)

exit
end if

;---------------------------------------------------------------------
; vertical integration: from surface to TOA
; u*D(h)/D(x)=[m/s]*[J/kg/m] = [W/kg]
; integral(u*D(h)/D(x))/g = [W/kg]*[kg/m/s2]/[m/s2] = [W/m2]
;---------------------------------------------------------------------
linlog = 1

; use Pa as the basic unit of pressure to integration: could lead to the result has unit of W/m2 by further multiplying "gravity-g".
p_bot2top			= lev(::-1)*100. ; hPa --> Pa
nlvl   				= dimsizes(p_bot2top) 
psfc   				= PS ; the unit of PS is Pa
pbot   				= 1100.*100.
ptop   				= 100.*100.
gravit				= 9.81

; caution: vibeta function needs at least three levels above surface.
; so doing land mask first could cause many warning information.
; solution: get OCNFRAC when reading 2D variables at the same time.

u_DhDx_1			= u_DhDx(case|:,lat|:,lon|:,lev|:)
u_DhDx_vint			= vibeta (p_bot2top,u_DhDx_1(:,:,:,::-1),linlog,psfc,pbot,ptop)
u_DhDx_vint			= u_DhDx_vint/gravit
copy_VarCoords_1(u_DhDx_1,u_DhDx_vint)

v_DhDy_1			= v_DhDy(case|:,lat|:,lon|:,lev|:)
v_DhDy_vint			= vibeta (p_bot2top,v_DhDy_1(:,:,:,::-1),linlog,psfc,pbot,ptop)
v_DhDy_vint			= v_DhDy_vint/gravit
copy_VarCoords_1(v_DhDy_1,v_DhDy_vint)

w_DhDp_1			= w_DhDp(case|:,lat|:,lon|:,lev|:)
w_DhDp_vint			= vibeta (p_bot2top,w_DhDp_1(:,:,:,::-1),linlog,psfc,pbot,ptop)
w_DhDp_vint			= w_DhDp_vint/gravit
copy_VarCoords_1(w_DhDp_1,w_DhDp_vint)

delete([/psfc/])
delete([/u_DhDx_1,v_DhDy_1,w_DhDp_1/])

;; NOTE: two questions left at 2017-11-15
; 1. DhDx has some problems: 
;    the extreme values mainly come from the large gradient over high topography, high latitudes and continent border.
;    Masking the land cover and restricting the region over tropics: 30S-30N, could solve this to some extent.
; 2. vertical integration result is not right.
;    the error comes from the unit of PS, which is Pa in default data.

;=====================================================================
; sum of all advection terms
;=====================================================================
advection			= u_DhDx_vint+v_DhDy_vint+w_DhDp_vint
copy_VarCoords(u_DhDx_vint,advection)

allsums			= advection + exforcing
copy_VarCoords(advection,allsums)

;=====================================================================
; write into NC file to check them
;=====================================================================
system("/bin/rm -f mse-budget-wo-obs_"+compare_type(icomp)+".nc")
ncdf = addfile("mse-budget-wo-obs_"+compare_type(icomp)+".nc","c")

ncdf->uwind			= uwind(:,:,{latS:latE},{lonS:lonE})
ncdf->vwind			= vwind(:,:,{latS:latE},{lonS:lonE})
ncdf->wwind			= wwind(:,:,{latS:latE},{lonS:lonE})
ncdf->MSE			= MSE(:,:,{latS:latE},{lonS:lonE})

ncdf->OCNFRAC			= OCNFRAC(:,{latS:latE},{lonS:lonE})

ncdf->DhDy			= DhDy(:,:,{latS:latE},{lonS:lonE})
ncdf->DhDx			= DhDx(:,:,{latS:latE},{lonS:lonE})
ncdf->DhDp			= DhDp(:,:,{latS:latE},{lonS:lonE})

ncdf->u_DhDx			= u_DhDx(:,:,{latS:latE},{lonS:lonE})
ncdf->v_DhDy			= v_DhDy(:,:,{latS:latE},{lonS:lonE})
ncdf->w_DhDp			= w_DhDp(:,:,{latS:latE},{lonS:lonE})

ncdf->u_DhDx_vint		= u_DhDx_vint(:,{latS:latE},{lonS:lonE})
ncdf->v_DhDy_vint		= v_DhDy_vint(:,{latS:latE},{lonS:lonE})
ncdf->w_DhDp_vint		= w_DhDp_vint(:,{latS:latE},{lonS:lonE})

ncdf->advection			= advection(:,{latS:latE},{lonS:lonE})

ncdf->LHFLX			= LHFLX(:,{latS:latE},{lonS:lonE})
ncdf->SHFLX			= SHFLX(:,{latS:latE},{lonS:lonE})
ncdf->RadCooling_allsky		= RadCooling_allsky(:,{latS:latE},{lonS:lonE})

ncdf->exforcing			= exforcing(:,{latS:latE},{lonS:lonE})

ncdf->allsums			= allsums(:,{latS:latE},{lonS:lonE})

delete([/uwind,vwind,wwind,MSE/])
delete([/OCNFRAC,PS/])
delete([/DhDy,DhDx,DhDp,u_DhDx,v_DhDy,w_DhDp,u_DhDx_vint,v_DhDy_vint,w_DhDp_vint,advection/])
delete([/LHFLX,SHFLX,RadCooling_allsky,exforcing,allsums/])

delete([/lat,lon,lev/])

end do ; do icomp=
exit
end if ; DataProcess

;=====================================================================
; Read data from mse-budget-wo-obs.nc and plotting
;=====================================================================

do icomp = 0,dimsizes(compare_type)-1

plot_vars			= (/"u_DhDx_vint","v_DhDy_vint","w_DhDp_vint","advection","exforcing","allsums"/)

ncin = addfile("mse-budget-wo-obs_"+compare_type(icomp)+".nc","r")
dSizes                          = getfiledimsizes(ncin)  
print(dSizes)

OCNFRAC				= ncin->OCNFRAC
printVarSummary(OCNFRAC)

tmp				= new((/dimsizes(plot_vars),dSizes(0),dSizes(2),dSizes(3)/),double) ; no lev dimension
tmp!0				= "vars"
tmp&vars			= ispan(0,dimsizes(plot_vars)-1,1)

delete([/dSizes/])

do ivar=0,dimsizes(plot_vars)-1
tmp(ivar,:,:,:)			= ncin->$plot_vars(ivar)$
end do

; mask land values
VarPlots			= tmp
VarPlots			= (/mask(tmp,OCNFRAC.ge.0.5,True)/)

printVarSummary(VarPlots)

; get difference between PDF and RH scheme
VarPlotsDiff			= VarPlots(:,1,:,:)-VarPlots(:,0,:,:)
copy_VarCoords_1(VarPlots(vars|:,lat|:,lon|:,case|:),VarPlotsDiff)
printVarSummary(VarPlotsDiff)

delete([/tmp,OCNFRAC/])

;------------------------------------------------------------------
; LonValuePlot
;------------------------------------------------------------------
latS				= -20.
latE				= -10.
lonS				= 120.
lonE				= 260.

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; get avarge over specific zonal bands
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lat_here 				= VarPlots&lat
rad 					= 4.0*atan(1.0)/180.
clat 					= cos(lat_here*rad)
clat!0 					= "lat"
clat&lat 				= lat_here

lon_here				= VarPlots&lon
lon_here!0				= "lon"
lon_here&lon				= VarPlots&lon
printVarSummary(lon_here)

VarPlotsZm				= dim_avg_wgt_n_Wrap(VarPlots(:,:,{latS:latE},:),clat({latS:latE}),1,2)
printVarSummary(VarPlotsZm)
VarPlotsDiffZm				= dim_avg_wgt_n_Wrap(VarPlotsDiff(:,{latS:latE},:),clat({latS:latE}),1,1)
printVarSummary(VarPlotsDiffZm)

if(LonValuePlot)then

wks				= gsn_open_wks(plot_type,"./figure/MSE-budget-wo-obs-LonValuePlot_"+compare_type(icomp)+"_"+latS+"_"+latE)

res				= True
res@gsnFrame			= False
res@gsnDraw			= False

res@vpWidthF			= 0.7
res@vpHeightF			= 0.25

res@tmXBLabelFontHeightF        = 0.02
res@tmYLLabelFontHeightF        = 0.02

res@xyLineThicknessF		= 3.0
res@xyLineColors		= (/"red","blue"/)

res@gsnYRefLine			= (/0.0,0.0/)
res@gsnYRefLineThicknessF	= 3.0
res@gsnYRefLineDashPattern	= 2

res@tiXAxisString		= "Longitude"
;res@tiYAxisString		= "Term"

res@gsnCenterStringFontHeightF	= 0.03

;res@gsnRightString		= "W/m~S~2~N~"
;res@gsnRightStringFontHeightF	= 0.03

xindexc				= ispan(tointeger(lonS),tointeger(lonE),20)
xlabelc				= where(xindexc.gt.180,xindexc-360.,xindexc)

res@tmXBMode			= "Explicit"
res@tmXBValues			= xindexc
res@tmXBLabels			= xlabelc

plot				= new(dimsizes(plot_vars),graphic)
plotDiff			= new(dimsizes(plot_vars),graphic)


res@trYMinF			= min(VarPlotsZm(:,:,{lonS:lonE}))
res@trYMaxF			= max(VarPlotsZm(:,:,{lonS:lonE}))
do ivar = 0,dimsizes(plot_vars)-1
	res@trYMinF			= min(VarPlotsZm(ivar,:,{lonS:lonE}))
	res@trYMaxF			= max(VarPlotsZm(ivar,:,{lonS:lonE}))
	res@gsnCenterString		= plot_vars(ivar)
	plot(ivar)			= gsn_csm_xy(wks,lon_here({lonS:lonE}),VarPlotsZm(ivar,:,{lonS:lonE}),res)
end do

resP				= True
resP@gsnPanelMainString		= "MSE budget from Park-RH"

gsn_panel(wks,plot,(/3,2/),resP)

;~~~~~~~~~~~~~~~~~~~~~
; Difference Plot
;~~~~~~~~~~~~~~~~~~~~~
res@trYMinF			= min(VarPlotsDiffZm(:,{lonS:lonE}))
res@trYMaxF			= max(VarPlotsDiffZm(:,{lonS:lonE}))

do ivar = 0,dimsizes(plot_vars)-1
	res@gsnCenterString		= plot_vars(ivar)
	plotDiff(ivar)			= gsn_csm_xy(wks,lon_here({lonS:lonE}),VarPlotsDiffZm(ivar,{lonS:lonE}),res)
end do

resP@gsnPanelMainString		= "MSE budget Difference b/t Gauss-PDF and Park-RH"

gsn_panel(wks,plotDiff,(/3,2/),resP)

delete(wks)
delete(plot)
delete(res)
delete(resP)


end if	; LonValuePlot

;------------------------------------------------------------------
; AreaMeanPlot
;------------------------------------------------------------------
latS				= -20.
latE				= -5.
lonS				= 240.
lonE				= 260.

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; get avarge over one specific region
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lat_here 				= VarPlots&lat
rad 					= 4.0*atan(1.0)/180.
clat 					= cos(lat_here*rad)
clat!0 					= "lat"
clat&lat 				= lat_here

lon_here				= VarPlots&lon
lon_here!0				= "lon"
lon_here&lon				= VarPlots&lon
printVarSummary(lon_here)

VarPlotsAm				= wgt_areaave_Wrap(VarPlots(:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}), 1.0, 0)
printVarSummary(VarPlotsAm)

VarPlotsDiffAm				= wgt_areaave_Wrap(VarPlotsDiff(:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,0)
printVarSummary(VarPlotsDiffAm)

if(AreaMeanPlot)then

wks				= gsn_open_wks(plot_type,"./figure/MSE-budget-wo-obs-AreaMeanPlot_"+compare_type(icomp)+"_"+latS+"_"+latE+"_"+lonS+"_"+lonE)

res				= True

res@gsnFrame			= False
res@gsnDraw			= False

res@vpWidthF			= 0.7
res@vpHeightF			= 0.25

;res@tmXBLabelFontHeightF        = 0.02
;res@tmYLLabelFontHeightF        = 0.02

colors				= (/"red","blue","brown"/)

res@xyMarkLineMode		= "Markers"
res@xyMarkers			= (/16,16/)
res@xyMarkerSizeF		= 0.015
res@xyMarkerColors		= colors

;res@tmXBLabelAngleF		= 45
;res@tmXBLabelFontHeightF	= 0.01

;res@tiXAxisString		= "Longitude"

;res@gsnRightString		= "W/m~S~2~N~"
;res@gsnStringFontHeightF	= 0.02

xindexc				= ispan(0,dimsizes(plot_vars)-1,1)
xlabelc				= (/"-uDhDx","-vDhDy","-wDhDp","advection","forcing","residual"/)

res@tmXBMode			= "Explicit"
res@tmXBValues			= xindexc
res@tmXBLabels			= xlabelc

miny				= min(VarPlotsAm)*1.1
maxy				= max(VarPlotsAm)*1.1
tmp				= new((/2/),double)
tmp				= (/abs(miny),abs(maxy)/)
xval				= max(tmp)
res@trYMinF			= -1.*xval
res@trYMaxF			= xval
delete(tmp)

plotAm				= gsn_csm_y(wks,transpose(VarPlotsAm),res)

;drawNDCGrid(wks)

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Difference plot
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
resDiff				= res
resDiff@xyMarkerColors		= (/"brown","brown","brown"/)

miny				= min(VarPlotsDiffAm)*1.1
maxy				= max(VarPlotsDiffAm)*1.1
tmp				= new((/2/),double)
tmp				= (/abs(miny),abs(maxy)/)
xval				= max(tmp)
resDiff@trYMinF			= -1.*xval
resDiff@trYMaxF			= xval
delete(tmp)

resDiff@tiYAxisString		= "Difference"

plotAmDiff			= gsn_csm_y(wks,VarPlotsDiffAm,resDiff)

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Draw all in the same plot
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plots				= gsn_csm_xy2(wks,xindexc,transpose(VarPlotsAm),VarPlotsDiffAm,res,resDiff)

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; add reference lines
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xpts				= (/0.0,dimsizes(plot_vars)-1/)
ypts				= (/0.0,0.0/)

res1				= True
res1@gsLineColor		= "grey"
res1@gsLineDashPattern		= 1
res1@gsLineThicknessF		= 3.0

dum1				= gsn_add_polyline(wks,plotAm,xpts,ypts,res1)
dum2				= gsn_add_polyline(wks,plotAmDiff,xpts,ypts,res1)
dum3				= gsn_add_polyline(wks,plots,xpts,ypts,res1)

;~~~~~~~~~~~~~~~~~~~~~
; add legends
;~~~~~~~~~~~~~~~~~~~~~
legends				= (/"Park-RH","Gauss-PDF","PDF-RH"/)

lgres				= True
lgres@lgMarkerColors		= colors
lgres@lgItemType		= res@xyMarkLineMode
lgres@lgMarkerIndexes		= (/16,16,16/)
lgres@lgLabelFontHeightF	= 0.08
lgres@vpWidthF			= 0.15
lgres@vpHeightF			= 0.1
lgres@lgPerimOn			= False

lbid1				= gsn_create_legend(wks,dimsizes(legends),legends,lgres)
lbid2				= gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres				= True
amres@amParallelPosF		= -0.40
amres@amOrthogonalPosF		= -0.30

annoid1				= gsn_add_annotation(plotAm,lbid1,amres)
annoid2				= gsn_add_annotation(plots,lbid2,amres)


gsn_panel(wks,(/plots/),(/1,1/),False)
;gsn_panel(wks,(/plotAm,plotAmDiff,plots/),(/3,1/),False)
gsn_panel(wks,(/plotAm,plotAmDiff/),(/2,1/),False)



delete(wks)
delete(plotAm)
delete(plotAmDiff)
delete(res)
delete(res1)
delete(resDiff)

end if	; LonValuePlot



;------------------------------------------------------------------
; LatLonPlot
;------------------------------------------------------------------
latS				= -30.
latE				= 30.
lonS				= 120.
lonE				= 290.


if(LatLonPlot)then

wks				= gsn_open_wks(plot_type,"./figure/MSE-budget-wo-obs-LatLonPlot_"+compare_type(icomp))

res				= True
res@gsnFrame			= False
res@gsnDraw			= False

res@cnFillOn			= True
res@cnLineLabelsOn		= False

res@mpMinLonF			= lonS
res@mpMaxLonF			= lonE
res@mpMinLatF			= latS
res@mpMaxLatF			= latE
res@mpCenterLonF	 	= 180.
res@gsnAddCyclic		= False
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Font height of tickmarks
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res@tmXBLabelFontHeightF	= 0.02
res@tmYLLabelFontHeightF	= 0.02
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Distance b/t labelbar and plot
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res@pmLabelBarOrthogonalPosF	= 0.20
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Labelbar
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;res@lbLabelBarOn		= False
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Font height of Strings
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res@gsnStringFontHeightF	= 0.03



plot				= new(dimsizes(plot_vars),graphic)
plotDiff			= new(dimsizes(plot_vars),graphic)

do ivar = 0,dimsizes(plot_vars)-1
	res@cnLevelSelectionMode	= "ManualLevels"
	res@cnMinLevelValF		= -200.
	res@cnMaxLevelValF		= 200.
	res@cnLevelSpacingF		= 20.

	res@gsnCenterString	= plot_vars(ivar)
	plot(ivar)		= gsn_csm_contour_map_ce(wks,VarPlots(ivar,0,:,:),res)
end do

resP				= True
resP@gsnPanelMainString		= "MSE budget from Park-RH"
;resP@gsnPanelLabelBar		= True

gsn_panel(wks,plot,(/3,2/),resP)
;~~~~~~~~~~~~~~~~~~~~~~~~~~
; Difference plots
;~~~~~~~~~~~~~~~~~~~~~~~~~~

do ivar = 0,dimsizes(plot_vars)-1
	res@cnLevelSelectionMode	= "ManualLevels"
	res@cnMinLevelValF		= -100.
	res@cnMaxLevelValF		= 100.
	res@cnLevelSpacingF		= 20.

	res@gsnCenterString	= plot_vars(ivar)
	plotDiff(ivar)		= gsn_csm_contour_map_ce(wks,VarPlotsDiff(ivar,:,:),res)
end do

resP@gsnPanelMainString		= "MSE budget Difference b/t Gauss-PDF and Park-RH"
gsn_panel(wks,plotDiff,(/3,2/),resP)

delete(wks)
delete(plot)
delete(res)
delete(resP)

end if	; LatLonPlot


delete([/clat,lat_here,lon_here/])
delete([/VarPlots,VarPlotsZm/])
delete([/VarPlotsDiff,VarPlotsDiffZm/])

end do; do icomp=
 
                              
end
