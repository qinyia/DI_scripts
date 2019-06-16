
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

function regrid_diff(var_obs[*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                     var_def[*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		     var_new[*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

;horizontal interpolation

lat_spc			= fspan(-90,90,72)
lon_spc			= fspan(0,360,144)

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

;var2D				= (/"SST","NetSfcHeat","FSNS","FLNS","LHFLX","SHFLX","CLDLOW"/)
;colors				= (/"black","brown","red","blue","green","orange","purple2"/)
;dashes				= fspan(0,dimsizes(var2D)-1,dimsizes(var2D))

;var2D				= (/"NetSfcHeat","FSNS","CLDLOW","CLDMED","CLDHGH","CLDTOT"/)
;colors				= (/"black","brown","burlywood3","bisque2","antiquewhite1"/)
;colors				= (/"black","brown","deepskyblue2","dodgerblue3","blue4","red"/)

;var2D				= (/"FSNS","PRECT","CLDHGH","CLDMED","CLDLOW","CLDTOT"/)
;colors				= (/"brown","blue","burlywood4","bisque3","antiquewhite2","darkgoldenrod4"/)
;dashes				= (/0,0,1,2,3,0/)

var2D				= (/"PRECT","CLDHGH","CLDMED","CLDLOW","CLDTOT"/)
colors				= (/"blue","burlywood4","bisque3","antiquewhite2","darkgoldenrod4"/)

dashes				= (/0,1,2,3,0/)

vector_L			= (/0/)
vector_R			= (/1,2,3,4/)

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
  var2D(ivar).eq."NetSfcHeat".or.var2D(ivar).eq."NetAtmHeat".or.var2D(ivar).eq."NetToaHeat")then
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

printVarSummary(var_interp)

;;;;; specific region mean; [0-20S, 180E-140W]

latS 					= -20.
latE 					= -5.
lonS					= 160.
lonE					= 280.


lat_here 				= var_interp&lat
rad 					= 4.0*atan(1.0)/180.
clat 					= cos(lat_here*rad)
clat!0 					= "lat"
clat&lat 				= lat_here

;========================================================
; Get meridional mean over (5S-10S) zonal band 
;========================================================
var_interp_mm		= dim_avg_wgt_n_Wrap(var_interp(:,{latS:latE},:),clat({latS:latE}),1,1)
var_interp_diff_mm	= dim_avg_wgt_n_Wrap(var_interp_diff(:,{latS:latE},:),clat({latS:latE}),1,1)

dims			= dimsizes(var_interp_mm)
print(dims)
if(ivar.eq.0)then
	var_interp_mm_all			= new((/dimsizes(var2D),dims(0),dims(1)/),float)
	var_interp_mm_all!0			= "vars"
	var_interp_mm_all(ivar,:,:)		= var_interp_mm
	printVarSummary(var_interp_mm_all)

	var_interp_diff_mm_all			= new((/dimsizes(var2D),dims(0),dims(1)/),float)
	var_interp_diff_mm_all!0		= "vars"
	var_interp_diff_mm_all(ivar,:,:)	= var_interp_diff_mm
	printVarSummary(var_interp_diff_mm_all)
else
	var_interp_mm_all(ivar,:,:)		= var_interp_mm
	var_interp_diff_mm_all(ivar,:,:)	= var_interp_diff_mm
end if

var_interp_mm_all&lon@units				= "degrees_east"
var_interp_diff_mm_all&lon@units			= "degrees_east"

;========================================================
; area mean over this specific region 
;========================================================

var_interp_areamean 			= wgt_areaave_Wrap(var_interp(:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1) ; obs, def, new
var_interp_diff_areamean 		= wgt_areaave_Wrap(var_interp_diff(:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1) ; def-obs,new-obs,new-def
print(var_interp_areamean)


if(ivar.eq.0)then
	var_interp_areamean_all 			= new((/dimsizes(var2D),3/),float)
	var_interp_areamean_all!0 			= "vars"
	var_interp_areamean_all!1 			= "cases"
	var_interp_areamean_all(ivar,:) 		= var_interp_areamean(:)

	var_interp_diff_areamean_all 			= new((/dimsizes(var2D),3/),float)
	var_interp_diff_areamean_all!0 			= "vars"
	var_interp_diff_areamean_all!1 			= "cases"
	var_interp_diff_areamean_all(ivar,:) 		= var_interp_diff_areamean(:)
else
	var_interp_areamean_all(ivar,:) 		= (/var_interp_areamean(:)/)
	var_interp_diff_areamean_all(ivar,:) 		= (/var_interp_diff_areamean(:)/)
end if


delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/var_interp,var_interp_diff/])
delete([/var_interp_mm,var_interp_diff_mm/])
delete(fapend)
delete(fdir)
delete([/lat_here,clat/])

end do ; do ivar=


;------------------------------------
;Plotting parameters
;------------------------------------

LonValuePlot			= True
AreaMeanPlot			= False

PlotType			= "pdf"

;-----------------------------------------------
; Plotting starts here...
;-----------------------------------------------

if(LonValuePlot)then

wks = gsn_open_wks(PlotType,"./figure/atm-LonValue-CLOUD_"+latS+"_"+latE)

plot				=new(2,graphic)

resL 				= True

resL@vpWidthF			= 0.7
resL@vpHeightF			= 0.4

resL@gsnFrame		 	= False
resL@gsnDraw			= False


resL@gsnYRefLine       		= 0.0
resL@gsnYRefLineColor		= "grey" 
resL@gsnYRefLineThicknessF	= 3.0
resL@gsnYRefLineDashPattern	= 3.0

resL@tiYAxisString		= "PRECT (mm/day)"

resL@tmXBTickSpacingF		= 15

resL@xyLineColors		= colors((/vector_L/))
resL@xyLineThicknessF 		= 4.0
resL@xyDashPatterns		= dashes((/vector_L/))

resL@trYMinF			= -30.
resL@trYMaxF			= -20.

resR				= True
resR@xyDashPatterns		= dashes((/vector_R/))
resR@xyLineColors		= colors((/vector_R/))
resR@xyLineThicknessF		= 4.0

resR@gsnYRefLine       		= 0.0
resR@gsnYRefLineColor		= "grey" 
resR@gsnYRefLineThicknessF	= 3.0
resR@gsnYRefLineDashPattern	= 3.0

resR@trYMinF			= -0.08
resR@trYMaxF			= -0.02

resR@tiYAxisString		= "Cloud Difference (fraction)"

tmp				= var_interp_mm_all(:,:,{lonS:lonE})

plot(0)			 	= gsn_csm_xy2(wks,tmp&lon,tmp((/vector_L/),0,:),tmp((/vector_R/),0,:),resL,resR) ; (vars,cases,lon)

delete(tmp)
;---------------------------------------------------------------
; Draw Lon-Value Plots: differences
;---------------------------------------------------------------

;resL@trYMinF			= -20.
;resL@trYMaxF			= 20.

tmp				= var_interp_diff_mm_all(:,:,{lonS:lonE}) ; vars,case,lon

plot(1)		 		= gsn_csm_xy2(wks,tmp&lon,tmp((/vector_L/),2,:),tmp((/vector_R/),2,:),resL,resR)

;--------------------------------------------
; start to add legends
;--------------------------------------------


legends                         = var2D
;legends				= (/"OBS","DEF","NEW"/)

lgres                           = True
lgres@lgLineColors            	= colors
lgres@lgLineThicknessF          = resL@xyLineThicknessF
lgres@lgDashIndexes		= dashes
lgres@lgLabelFontHeightF        = 0.08
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.15
lgres@lgPerimOn                 = False
lgres@lgItemOrder		= fspan(dimsizes(legends)-1,0,dimsizes(legends))

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)
lbid2                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres                           = True
;amres@amJust			= "BottomLeft"
;amres@amParallelPosF            = -0.50
;amres@amOrthogonalPosF          = 0.50

amres@amJust			= "TopLeft"
amres@amParallelPosF            = -0.5
amres@amOrthogonalPosF          = -0.5

annoid1                         = gsn_add_annotation(plot(0),lbid1,amres)
annoid2                         = gsn_add_annotation(plot(1),lbid2,amres)

gsn_panel(wks,plot(1),(/1,1/),False)

;gsn_panel(wks,plot,(/2,1/),False)

delete(resL)
delete(plot)
delete(wks)

end if ; LonValuePlot


if(AreaMeanPlot)then
;==============================================================================================================
wks 					= gsn_open_wks(PlotType,"./figure/atm_energy_budget_"+compare_type(icomp)+"_"+latS+"_"+latE+"_"+lonS+"_"+lonE)

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
cmap2 					= read_colormap_file("ViBlGrWhYeOrRe")

res                       		= True                 ; plot mods desired         
res@gsnFrame              		= False                ; don't advance frame yet

res@tmXMajorGrid 			= True
res@tmXMajorGridThicknessF 		= 1.0
res@tmXMajorGridLineDashPattern		= 2

res@tmYMajorGrid 			= True
res@tmYMajorGridThicknessF 		= 1.0
res@tmYMajorGridLineDashPattern		= 2

res@gsnYRefLine 			= 0.0

res@tmXBLabelFontHeightF		= 0.03
res@tmYLLabelFontHeightF		= 0.02


res@gsnXYBarChart         		= True                 ; turn on bar chart
res@gsnXYBarChartBarWidth 		= 2.0

;res@tmXBOn                		= False                ; turn off tickmarks at bot

res@trYMinF 				= min(var_interp_areamean_all)-5.0
res@trYMaxF 				= max(var_interp_areamean_all)+5.0
res@trXMinF 				= 0.
res@trXMaxF 				= 50.0

res@tiMainString          		= "Heat budget"


;x1 					= fspan(2.0,40.0,dimsizes(var2D))
;x2 					= fspan(3.0,41.0,dimsizes(var2D))
;x3 					= fspan(4.0,42.0,dimsizes(var2D))

x1 					= fspan(2.0,40.0,dimsizes(var2D))
x2 					= fspan(4.0,42.0,dimsizes(var2D))
x3 					= fspan(6.0,44.0,dimsizes(var2D))



res@tmXBMode   				= "Explicit"
res@tmXBValues 				= x2
res@tmXBLabels 				= var2D

res@tmXBLabelAngleF 			= 45
res@tmXBLabelFontHeightF 		= 0.01


res@gsnXYBarChartColors2 		= (/"grey"/)
plot1 					= gsn_csm_xy(wks,x1,var_interp_areamean_all(:,0),res)
res@gsnXYBarChartColors2 		= (/"red"/)
plot2 					= gsn_csm_xy(wks,x2,var_interp_areamean_all(:,1),res)
res@gsnXYBarChartColors2 		= (/"blue"/)
plot3 					= gsn_csm_xy(wks,x3,var_interp_areamean_all(:,2),res)

;**********************************************************
; add labelbar to second plot
;********************************************************** 
  getvalues plot1                ; get plot size for use in creating labelbar
  "vpXF"      : vpx
  "vpYF"      : vpy
  "vpHeightF" : vph
  "vpWidthF"  : vpw
  end getvalues

  lbw    = 0.2 * vpw           ; Make labelbar size a fraction of the plot.
  lbh    = 0.25 * vph
  nboxes = 3

labels 					= (/"OBS","DEF","NEW"/)
   
  lbres                    		= True          ; labelbar only resources
  lbres@vpWidthF           		= 0.2 * vpw     ; labelbar width
  lbres@vpHeightF          		= 0.25 * vph     ; labelbar height
  lbres@lbBoxMajorExtentF  		= 0.75          ; puts space between color boxes
  lbres@lbFillColors       		= (/"grey","red","blue"/) ; labelbar colors
  lbres@lbMonoFillPattern  		= True          ; Solid fill pattern
  lbres@lbLabelFontHeightF 		= 0.008         ; font height. default is small
  lbres@lbLabelJust        		= "CenterLeft"  ; left justify labels
   
  gsn_labelbar_ndc(wks,nboxes,labels,vpx+vph-lbw,0.77,lbres)
   
draw(plot1) 
frame(wks)


res@trYMinF 				= min(var_interp_diff_areamean_all)-5.0 
res@trYMaxF 				= max(var_interp_diff_areamean_all)+5.0


res@gsnXYBarChartColors2 		= (/"grey"/)
plot3 					= gsn_csm_xy(wks,x1,var_interp_diff_areamean_all(:,0),res)

res@gsnXYBarChartColors2 		= (/"red"/)
plot4 					= gsn_csm_xy(wks,x2,var_interp_diff_areamean_all(:,1),res)

res@gsnXYBarChartColors2 		= (/"blue"/)
plot5 					= gsn_csm_xy(wks,x3,var_interp_diff_areamean_all(:,2),res)

labels 					= (/"DEF-OBS","NEW-OBS","NEW-DEF"/)
gsn_labelbar_ndc(wks,nboxes,labels,vpx+vph-lbw,0.77,lbres)

draw(plot3)
frame(wks)

end if ; AreaMeanPlot


end do; do icomp=
 
                              
end
