
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

;var2D				= (/"SST","NetSfcHeat","FSNS","FLNS","LHFLX","SHFLX","CLDLOW","LHFLX_revised"/)
;colors				= (/"black","brown","red","blue","green","orange","purple","gray"/)
;thick				= (/5,5,3,3,3,3,5,3/)
;dashes				= (/0,0,14,14,14,14,0,14/)

var2D				= (/"SST","NetSfcHeat","FSNS","FLNS","LHFLX","SHFLX","CLDLOW"/)
colors				= (/"black","brown","red","blue","green","orange","purple"/)
thick				= (/5,5,3,3,3,3,5/)
dashes				= (/0,0,14,14,14,14,0/)

;var2D				= (/"SST","NetSfcHeat","FSNS","FLNS","LHFLX_revised","SHFLX","CLDLOW"/)
;colors				= (/"black","brown","red","blue","green","orange","purple"/)
;thick				= (/5,5,3,3,3,3,5/)
;dashes				= (/0,0,14,14,14,14,0/)

;var2D				= (/"NetSfcHeat","NetAtmHeat","NetToaHeat"/)
;colors				= (/"black","brown","red"/)


;compare_type 			= (/"amip-obs"/)
compare_type 			= (/"cmip-obs"/)


case				= (/"def","new"/)

endyear				= (/"0004","0006","0008","0012"/)

Dir				= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

;sstvt_case		= (/"BC5_f19g16_mac2_cosp","BC5_f19g16_sstvt_0.52",\
;			    "BC5_f19g16_sstvt_enlarge_region_0.75","BC5_f19g16_sstvt_abs_100percent","BC5_f19g16_sstvt_0.85"/)
sstvt_case		= (/"BC5_f19g16_sstvt_0.85"/)
;sstvt_case		= (/"BC5_f19g16_sstvt_enlarge_region_0.75"/)
;sstvt_case		= (/"BC5_f19g16_sstvt_abs_100percent"/)

do isstcase = 0,dimsizes(sstvt_case)-1

do iend = 0,dimsizes(endyear)-1

do ivar=0,dimsizes(var2D)-1
print(var2D(ivar))

if(compare_type(icomp) .eq."cmip-obs")then
	casename = (/"","BC5_f19g16_cosp",sstvt_case(isstcase)/)
	fdir = (/Dir+"obs/",Dir+"atm/",Dir+"atm/"/)

end if

fapend 					= (/"_ANN_climo.nc", \
       					    "_ANN_climo_"+endyear(iend)+".nc", \
	   				    "_ANN_climo_"+endyear(iend)+".nc" /)


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
  var2D(ivar).eq."LHFLX_revised")then
	casename(0) 			= "JRA25"
end if


print(fdir(0))
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
lonS					= 178.5
lonE					= 262.5

;latS 					= -90.
;latE 					= 90.
;lonS					= 0.
;lonE					= 360.


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


; according to the comment from Reviewer 2, double check the the net sfc heat flux and the sum of its components.
if(var2D(dimsizes(var2D)-1).eq."LHFLX_revised")then

Net_mm						= var_interp_mm_all(1,:,:)
FSNS_mm						= var_interp_mm_all(2,:,:)
FLNS_mm						= var_interp_mm_all(3,:,:)
LH_mm						= var_interp_mm_all(4,:,:)
SH_mm						= var_interp_mm_all(5,:,:)
LHR_mm						= var_interp_mm_all(7,:,:)

Sum1						= FSNS_mm-FLNS_mm-LH_mm-SH_mm
Sum2						= FSNS_mm-FLNS_mm-LHR_mm-SH_mm
copy_VarCoords(FSNS_mm,Sum1)
copy_VarCoords(FSNS_mm,Sum2)

lonsp						= 165.

; write table.
header						= (/"--------------------------------------------------------",\
					            "This is an output to verify sfc energy balance at "+lonsp,\
			     			    "--------------------------------------------------------"/)
hlist						= [/header/]

fname						= "fig4-sfc-energy-check-at-"+lonsp+".txt"

a						= (/"JRA25","Park-RH","Gauss-PDF"/)
b						= (/Sum1(0,{lonsp}),Sum1(1,{lonsp}),Sum1(2,{lonsp})/)
c						= (/Sum2(0,{lonsp}),Sum2(1,{lonsp}),Sum2(2,{lonsp})/)
d						= (/Net_mm(0,{lonsp}),Net_mm(1,{lonsp}),Net_mm(2,{lonsp})/)
write_table(fname,"w",hlist,"%s")
write_table(fname,"a",[/"SUM1-LH_O","SUM2-LH_N","NetSfcHeat"/],"	%s 	%s 	%s")
write_table(fname,"a",[/a,b,c,d/],"%s %9.3f %9.3f %9.3f")

end if


;------------------------------------
;Plotting parameters
;------------------------------------

LonValuePlot			= True
AreaMeanPlot			= False

PlotType			= "x11"

;-----------------------------------------------
; Plotting starts here...
;-----------------------------------------------

if(LonValuePlot)then
if(var2D(dimsizes(var2D)-1).eq."LHFLX_revised")then
	wks = gsn_open_wks(PlotType,"./figure/fig4-atm-heat-budget-LonValue_"+latS+"_"+latE+"_refined_2LHFLX")
else
	wks = gsn_open_wks(PlotType,"./figure/fig4-atm-heat-budget-LonValue_"+latS+"_"+latE+"_refined_revision_1_"+sstvt_case(isstcase)+"_"+endyear(iend))
end if

plot				=new(2,graphic)

resL 				= True

resL@vpWidthF			= 0.7
resL@vpHeightF			= 0.45

resL@gsnFrame		 	= False
resL@gsnDraw			= False


resL@gsnYRefLine       		= 0.0
resL@gsnYRefLineColor		= "black" 
resL@gsnYRefLineThicknessF	= 2.0
resL@gsnYRefLineDashPattern	= 0.0

;resL@tiYAxisString		= "heat flux [W/m~S~2~N~] / cloud [%] DIFF" 
resL@tiYAxisString		= "heat flux (W/m~S~2~N~) / cloud (%) DIFF" 

;resL@tmXBMode        		= "Manual"	
;resL@tmXBTickStartF  		= lonS
;resL@tmXBTickEndF    		= lonE
resL@tmXBTickSpacingF		= 15

;resL@tiXAxisString		= "longitude (Degree_east)"

resL@xyLineColors		= colors(1:)
resL@xyLineThicknesses 		= thick(1:)
resL@xyDashPatterns		= dashes(1:)

resL@trYMinF			= -30.
resL@trYMaxF			= 30.


resR				= True
resR@xyDashPatterns		= dashes(0)
resR@xyLineColors		= colors(0)
resR@xyLineThicknessF		= thick(0)

resR@trYMinF			= -2.
resR@trYMaxF			= 2.

resR@tiYAxisString		= "SST DIFF (~S~o~N~C)"

tmp				= var_interp_mm_all(:,:,{lonS:lonE})

plot(0)			 	= gsn_csm_xy2(wks,tmp&lon,tmp(1:,0,:),tmp(0,0,:),resL,resR) ; (vars,cases,lon)

delete(tmp)
;---------------------------------------------------------------
; Draw Lon-Value Plots: differences
;---------------------------------------------------------------

;resL@trYMinF			= -20.
;resL@trYMaxF			= 20.

tmp				= var_interp_diff_mm_all(:,:,{lonS:lonE})

plot(1)		 		= gsn_csm_xy2(wks,tmp&lon,tmp(1:,2,:),tmp(0,2,:),resL,resR)

;--------------------------------------------
; start to add legends
;--------------------------------------------


if(.not.any(var2D.eq."LHFLX").and.\ ; no LHFLX in var2D
    any(var2D.eq."LHFLX_revised").and.\ ; yes LHFLX_revised in var2D
    var2D(dimsizes(var2D)-1).ne."LHFLX_revised")\ ; LHFLX_revised is not the last one
    then
	legends                         = where(var2D.eq."LHFLX_revised","LHFLX",var2D)
else
	legends				= var2D

end if

;legends				= (/"OBS","DEF","NEW"/)

lgres                           = True
lgres@lgLineColors            	= colors
lgres@lgLineThicknessF          = 3.0
lgres@lgDashIndexes		= dashes
lgres@lgLabelFontHeightF        = 0.08
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.15
lgres@lgPerimOn                 = False
lgres@lgItemOrder		= fspan(dimsizes(legends)-1,0,dimsizes(legends))

lbid1                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)
lbid2                           = gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres                           = True
amres@amJust			= "BottomLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = 0.50

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

res@gsnYRefLine 			= 0.0

res@tmXBLabelFontHeightF		= 0.02
res@tmYLLabelFontHeightF		= 0.02
res@tmYRLabelFontHeightF		= 0.02


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

end do ; do iend=
end do ; do isstcase=


end do; do icomp=
 
                              
end
