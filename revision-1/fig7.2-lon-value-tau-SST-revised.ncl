
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

function regrid_diff(var_obs[*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                     var_def[*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		     var_new[*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

;horizontal interpolation

lat_spc			= fspan(-90,90,72)
lon_spc			= fspan(0,360,181)

var_interp 				= lat_lon_interpolate_3vars_specgrid(var_obs,lat_obs,lon_obs, \
				       				             var_def,lat_def,lon_def, \
				       				             var_new,lat_new,lon_new,\
									     lat_spc,lon_spc)
var_interp_diff 			= var_interp

var_interp_diff 			= var_interp
var_interp_diff(0,:,:) 			= var_interp(1,:,:) - var_interp(0,:,:) ; def-obs
var_interp_diff(1,:,:) 			= var_interp(2,:,:) - var_interp(0,:,:) ; new-obs
var_interp_diff(2,:,:) 			= var_interp(2,:,:) - var_interp(1,:,:) ; new-def

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin

;var2D				= (/"FSNS","PRECT","CLDHGH","CLDMED","CLDLOW","CLDTOT"/)
;colors				= (/"brown","blue","burlywood4","bisque3","antiquewhite2","darkgoldenrod4"/)
;dashes				= (/0,0,1,2,3,0/)

var2D				= (/"SST","TAUX"/)
colors				= (/(/"blue","blue"/),(/"burlywood4","burlywood4"/)/)
;colors				= (/(/"blue","blue"/),(/"red","red"/)/)

colors_1d			= (/"blue","burlywood4"/)

dashes				= (/(/0,1/),(/0,1/)/)

;dashes				= fspan(0,dimsizes(var2D)-1,dimsizes(var2D))


;compare_type 			= (/"amip-obs"/)
compare_type 			= (/"cmip-obs"/)


case				= (/"def","new"/)

Dir				= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

do ivar=0,dimsizes(var2D)-1
print(var2D(ivar))


	if(compare_type(icomp) .eq."amip-obs")then
		if(var2D(ivar).eq."SST")then
			print("cannot process SST in this compare_type:"+compare_type(icomp))
			print("please remove SST in varible list")
			exit
		end if
		casename = (/"","FC5_f09f09_MG15_COSP","new-cldfrac16-1d-MG15-revise-2-COSP"/)
		fdir = (/Dir+"obs/",Dir+"amip/",Dir+"amip/"/)

	end if
	if(compare_type(icomp) .eq."cmip-obs")then
		casename = (/"","BC5_f19g16","BC5_f19g16_mac2"/)
		fdir = (/Dir+"obs/",Dir+"atm/",Dir+"atm/"/)

	end if


fapend 					= (/"_ANN_climo.nc", \
       					    "_ANN_climo.nc", \
			                    "_ANN_climo.nc" /)


if(var2D(ivar).eq."SST")then
	casename(0) 			= "HadISST_PD"
	fdir 				= (/Dir+"obs/", \
	     				    Dir+"ocn/", \
		 			    Dir+"ocn/"/)
end if

if(var2D(ivar).eq."PRECT")then
	casename(0) 			= "GPCP"
end if

if(var2D(ivar).eq."CLDLOW".or.var2D(ivar).eq."CLDMED".or.var2D(ivar).eq."CLDHGH".or.var2D(ivar).eq."CLDTOT")then
	casename(0) 			= "CLOUDSAT"
end if

if(var2D(ivar).eq."FSNS" .or. var2D(ivar).eq."FLNS" .or. var2D(ivar).eq."LHFLX" .or. var2D(ivar).eq."SHFLX".or.\
  var2D(ivar).eq."FSNT".or.var2D(ivar).eq."FLUT".or.\
  var2D(ivar).eq."NetSfcHeat".or.var2D(ivar).eq."NetAtmHeat".or.var2D(ivar).eq."NetToaHeat".or.\
  var2D(ivar).eq."TAUX".or.var2D(ivar).eq."TAUY")then
	casename(0) 			= "JRA25"
end if



result_obs				= get_var_2D(fdir(0),fapend(0),casename(0),var2D(ivar),"obs")
var_obs					= result_obs[0]
lat_obs					= result_obs[1]
lon_obs					= result_obs[2]

result_def				= get_var_2D(fdir(1),fapend(1),casename(1),var2D(ivar),"model")
var_def					= result_def[0]
lat_def					= result_def[1]
lon_def					= result_def[2]

result_new				= get_var_2D(fdir(2),fapend(2),casename(2),var2D(ivar),"model")
var_new					= result_new[0]
lat_new					= result_new[1]
lon_new					= result_new[2]

; note 1: regridding will not affect the global/regional mean values significantly.
; note 2: have verified that global mean atmosphere budget (cmip and amip): it is conserved in OBS and two simulations.

result 					= regrid_diff(var_obs,lat_obs,lon_obs,\
       					              var_def,lat_def,lon_def,\
		              		              var_new,lat_new,lon_new)

var_interp 				= result[0]
var_interp_diff 			= result[1]

; change units
if(var2D(ivar).eq."TAUX".or.var2D(ivar).eq."TAUY")then
	var_interp			= var_interp*10.	; N/m2 -> dyn/cm2
	var_interp_diff			= var_interp_diff*10.
end if

printVarSummary(var_interp)

;========================================================
; Save all variables into one array
;========================================================
dims						= dimsizes(var_interp)
if(ivar.eq.0)then
	vars_2D				= new((/dimsizes(var2D),dims(0),dims(1),dims(2)/),float)
	vars_2D!0			= "vars"
	vars_2D(ivar,:,:,:)		= var_interp

	vars_2D_diff			= new((/dimsizes(var2D),dims(0),dims(1),dims(2)/),float)
	vars_2D_diff!0			= "vars"
	vars_2D_diff(ivar,:,:,:)	= var_interp_diff
else
	vars_2D(ivar,:,:,:)		= var_interp
	vars_2D_diff(ivar,:,:,:)		= var_interp_diff
end if

vars_2D&lon@units			= "degrees_east"
vars_2D_diff&lon@units			= "degrees_east"

printVarSummary(vars_2D)

delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/var_interp,var_interp_diff/])
delete(fapend)
delete(fdir)
end do ; do ivar=

;========================================================
; calculate horizontal advection: zonal and meridional
;========================================================
; get latitudinal and longitudinal gradient
SST	 		= 	vars_2D(0,:,:,:)
printVarSummary(SST)

gradLatLon  		= 	grad_latlon_cfd (SST, vars_2D&lat, vars_2D&lon, True, False)
                	               
T_grad_dy  		= 	gradLatLon[0] 
T_grad_dx  		= 	gradLatLon[1]

delete(gradLatLon)

; change to 10^-6 K/m
T_grad_dx		= 	T_grad_dx*1e6

T_grad_dx@long_name 	= 	"zonal gradient"
T_grad_dy@long_name 	= 	"meridional gradient"
T_grad_dx@units     	= 	"degC/m"
T_grad_dy@units     	= 	"degC/m"

printVarSummary(T_grad_dy)

;====================================================================
; the difference of T_grad_dx
;====================================================================
T_grad_dx_diff		= T_grad_dx
T_grad_dx_diff(0,:,:)	= T_grad_dx(1,:,:)-T_grad_dx(0,:,:)
T_grad_dx_diff(1,:,:)	= T_grad_dx(2,:,:)-T_grad_dx(0,:,:)
T_grad_dx_diff(2,:,:)	= T_grad_dx(2,:,:)-T_grad_dx(1,:,:)

printVarSummary(T_grad_dx_diff)
printMinMax(T_grad_dx_diff,True)

; add T_grad_dx into the whole big array
vars_2D_new					= new((/dimsizes(var2D)+1,dims(0),dims(1),dims(2)/),float)
vars_2D_new!0					= "vars"
vars_2D_new(:dimsizes(var2D)-1,:,:,:)		= vars_2D
vars_2D_new(dimsizes(var2D),:,:,:)		= T_grad_dx
printVarSummary(vars_2D_new)

vars_2D_diff_new					= new((/dimsizes(var2D)+1,dims(0),dims(1),dims(2)/),float)
vars_2D_diff_new!0					= "vars"
vars_2D_diff_new(:dimsizes(var2D)-1,:,:,:)		= vars_2D_diff
vars_2D_diff_new(dimsizes(var2D),:,:,:)			= T_grad_dx_diff
printVarSummary(vars_2D_diff_new)

;;;;; specific region mean; [0-20S, 180E-140W]

latS 					= -20.
latE 					= -5.
lonS					= 180.
lonE					= 260.


lat_here 				= vars_2D&lat
rad 					= 4.0*atan(1.0)/180.
clat 					= cos(lat_here*rad)
clat!0 					= "lat"
clat&lat 				= lat_here

;========================================================
; Get meridional mean over (5S-10S) zonal band 
;========================================================
; vars,case,lat,lon
vars_2D_mm		= dim_avg_wgt_n_Wrap(vars_2D_new(:,:,{latS:latE},:),clat({latS:latE}),1,2)
vars_2D_diff_mm	= dim_avg_wgt_n_Wrap(vars_2D_diff_new(:,:,{latS:latE},:),clat({latS:latE}),1,2)
printVarSummary(vars_2D_mm)

;========================================================
; area mean over this specific region 
;========================================================

vars_2D_areamean 			= wgt_areaave_Wrap(vars_2D_new(:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1) ; obs, def, new
vars_2D_diff_areamean 		= wgt_areaave_Wrap(vars_2D_diff_new(:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1) ; def-obs,new-obs,new-def
printVarSummary(vars_2D_areamean)

vector_L				= (/0/)
vector_R				= (/1/)

var2D_new				= (/"SST","TAUX","dSSTdx"/)
colors_new				= (/(/"blue","blue"/),(/"burlywood4","burlywood4"/),(/"red","red"/)/)
colors_1d_new				= (/"blue","burlywood4","red"/)
thick_new				= (/(/2.0,2.0/),(/2.0,2.0/),(/2.0,2.0/)/)
dashes_new				= (/(/0,14/),(/0,14/),(/0,14/)/)

;------------------------------------
;Plotting parameters
;------------------------------------

LonValuePlot			= True
AreaMeanPlot			= False

;-----------------------------------------------
; Plotting starts here...
;-----------------------------------------------

if(LonValuePlot)then

wks = gsn_open_wks("pdf","./figure/fig7.2-atm-LonValue-SST-tau_"+latS+"_"+latE+"_revised_revision_1")

plot				=new(3,graphic)

res 				= True
res@vpWidthF			= 0.7
res@vpHeightF			= 0.4
res@gsnFrame		 	= False
res@gsnDraw			= False

res@tmXBTickSpacingF		= 15
res@tmXBLabelFontHeightF	= 0.0325
res@tmYLLabelFontHeightF	= 0.0325
res@tmYRLabelFontHeightF	= 0.0325

;res@tmLabelAutoStride		= True
res@tmYLLabelStride		= 2
res@tmYRLabelStride		= 2

;res@tmXBMode			= "Explicit"
;res@tmXBValues			= (/160,180,200,220,240,260,280/)
;res@tmXBLabels			= (/"160E","180","160W","140W","120W","100W","80W"/)

; plot zonal SSTand zonal wind stress

resL 				= res

; warning: these two lines will give out some warning messages about get_color_name.
; I don't know the reason for this now.
resL@gsnYRefLine       		= 0.0
resL@gsnYRefLineColor		= "grey"

resL@xyLineColors		= colors_new((/vector_L/),:)
resL@xyLineThicknesses 		= thick_new((/vector_L/),:)
resL@xyDashPatterns		= dashes_new((/vector_L/),:)
resL@tmYLLabelFontColor 	= resL@xyLineColors

resL@trYMinF			= 23.
resL@trYMaxF			= 30.

resR				= res

resR@xyLineColors		= colors_new((/vector_R/),:)
resR@xyDashPatterns		= dashes_new((/vector_R/),:)
resR@xyLineThicknesses		= thick_new((/vector_R/),:)
resR@tmYRLabelFontColor 	= resR@xyLineColors

resR@gsnYRefLine       		= 0.0
resR@gsnYRefLineColor		= "grey" 

resR@trYMinF			= -0.8
resR@trYMaxF			= -0.2


;resL@tiYAxisString		= var2D_new((/vector_L/))
;resL@tiYAxisString		= "SST"
resL@tiYAxisString		= ""
resR@tiYAxisString		= ""

tmp				= vars_2D_mm(:,:,{lonS:lonE})

; obs and def (/0,1/)
plot(0)			 	= gsn_csm_xy2(wks,tmp&lon,tmp((/vector_L/),(/0,1/),:),tmp((/vector_R/),(/0,1/),:),resL,resR) ; (vars,case,lon)

resL@tiYAxisString		= ""
resR@tiYAxisString		= ""
; obs and new (/0,2/)
plot(1)			 	= gsn_csm_xy2(wks,tmp&lon,tmp((/vector_L/),(/0,2/),:),tmp((/vector_R/),(/0,2/),:),resL,resR) ; (vars,case,lon)

delete(tmp)
;---------------------------------------------------------------
; Draw Lon-Value Plots: differences
;---------------------------------------------------------------


resL@xyLineThicknesses 		= thick_new((/vector_L,vector_R/),(/1/))
resR@xyLineThicknesses 		= thick_new((/vector_L,vector_R/),(/1/))
resL@xyDashPatterns		= dashes_new((/vector_L,vector_R/),(/1/))
resR@xyDashPatterns		= dashes_new((/vector_L,vector_R/),(/1/))

resL@tiYAxisString		= ""
;resR@tiYAxisString		= var2D_new((/vector_R/))
resR@tiYAxisString		= ""
;resR@tiYAxisString		= "TAUX"
;resR@tiYAxisFontColor		= resR@xyLineColors
;resR@tiYAxisFontHeightF		= 0.02

resL@trYMinF			= -2.0
resL@trYMaxF			= 0.

resR@trYMinF			= -0.24
resR@trYMaxF			= 0.


tmp				= vars_2D_diff_mm(:,:,{lonS:lonE}) ; vars,case,lon

plot(2)		 		= gsn_csm_xy2(wks,tmp&lon,tmp((/vector_L/),2,:),tmp((/vector_R/),2,:),resL,resR)

;--------------------------------------------
; start to add legends
;--------------------------------------------


legends                         = var2D_new(:1)
;legends				= (/"OBS","DEF","NEW"/)

lgres                           = True
lgres@lgLineColors            	= colors_1d_new((/vector_L,vector_R/))
lgres@lgLineThicknessF          = 3.0
lgres@lgDashIndexes		= (/0,0/)
lgres@lgLabelFontHeightF        = 0.20
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.20
lgres@lgPerimOn                 = False
lgres@lgItemOrder		= fspan(dimsizes(legends)-1,0,dimsizes(legends))

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)
lbid2                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres                           = True
amres@amJust			= "BottomLeft"
amres@amParallelPosF            = -0.30
amres@amOrthogonalPosF          = 0.50

;amres@amJust			= "TopLeft"
;amres@amParallelPosF            = -0.5
;amres@amOrthogonalPosF          = -0.5

annoid1                         = gsn_add_annotation(plot(0),lbid1,amres)
;annoid2                         = gsn_add_annotation(plot(1),lbid2,amres)

resP				= True
resP@gsnPanelXWhiteSpacePercent	= 0
resP@gsnPanelYWhiteSpacePercent	= 0

resP@amJust			= "BottomLeft"
resP@gsnPanelFigureStrings      = (/"g)","h)","i)"/)

;gsn_panel(wks,plot,(/3,1/),False)
gsn_panel(wks,plot,(/1,3/),resP)


delete(resL)
delete(plot)
delete(wks)
exit

end if ; LonValuePlot


end do; do icomp=
 
                              
end
