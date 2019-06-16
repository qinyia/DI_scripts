
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

;var2D				= (/"CrfToa","CrfSfc","Acrf"/)
;var2D				= (/"SWCF_Toa","LWCF_Toa","SWCF_Sfc","LWCF_Sfc","CrfToa","CrfSfc"/)
;var2D				= (/"CrfToa","CrfSfc","ACRF"/)
var2D				= (/"RadCooling_allsky","RadCooling_clrsky","ACRF"/)

var3D				= (/"MSE"/)


;compare_type 			= (/"cmip-obs"/)
compare_type 			= (/"amip-obs"/)


Dir				= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

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

	printVarSummary(var_def)
	printMinMax(var_def,True)

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
	
end do ; do ivar=

if(False)then
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

	result 					= regrid_diff_3D(var_def,lev_def,lat_def,lon_def,\
		              		                         var_new,lat_new,lon_new)

	var_3D 				= result[0]
	var_3D_diff 			= result[1]

	delete([/result,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])

	; put all 2D variables into one variable

	index				= dimsizes(var_3D)
	if(ivar.eq.0)then
		vars_LevLatLon		= new((/dimsizes(var3D),index(0),index(1),index(2),index(3)/),float)
		vars_LevLatLon!0		= "vars"
		vars_LevLatLon&vars	= ispan(0,dimsizes(var3D)-1,1)
	end if

	delete(index)
	vars_LevLatLon(ivar,:,:,:,:)		= var_3D

end do; do ivar=

printVarSummary(vars_LevLatLon)
printVarSummary(vars_LatLon)

end if

;====================================================
; initially test  the basic feature of these variables
;====================================================

if(True)then
wks			= gsn_open_wks("x11","./figure/"+var2D(0))

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




;;;;; specific region mean; [0-20S, 180E-140W]

latS 					= -20.
latE 					= 0.
lonS					= 180.
lonE					= 220.


lat_here 				= var_interp&lat
rad 					= 4.0*atan(1.0)/180.
clat 					= cos(lat_here*rad)
clat!0 					= "lat"
clat&lat 				= lat_here

; area mean over this specific region 
var_interp_areamean 			= wgt_areaave_Wrap(var_interp(:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1) ; obs, def, new
var_interp_diff_areamean 		= wgt_areaave_Wrap(var_interp_diff(:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1) ; def-obs,new-obs,new-def

print(var_interp_areamean)

delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/var_interp,var_interp_diff/])
delete(fapend)
delete(fdir)
delete([/lat_here,clat/])

if(ivar.eq.0)then
	var_interp_all 			= new((/dimsizes(var2D),3/),float)
	var_interp_all!0 		= "vars"
	var_interp_all&vars 		= var2D
	var_interp_all!1 		= "cases"
	var_interp_all&cases 		= (/"obs","def","new"/)
	var_interp_all(ivar,:) 		= (/var_interp_areamean(:)/)

	var_interp_diff_all 		= new((/dimsizes(var2D),3/),float)
	var_interp_diff_all!0 		= "vars"
	var_interp_diff_all&vars 	= var2D
	var_interp_diff_all!1 		= "cases"
	var_interp_diff_all&cases 	= (/"obs","def","new"/)
	var_interp_diff_all(ivar,:) 	= (/var_interp_diff_areamean(:)/)
else
	var_interp_all(ivar,:) 		= (/var_interp_areamean(:)/)
	var_interp_diff_all(ivar,:) 	= (/var_interp_diff_areamean(:)/)
end if



;==============================================================================================================
wks 					= gsn_open_wks("x11","./figure/surface_energy_budget_"+compare_type(icomp))

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


res@gsnXYBarChart         		= True                 ; turn on bar chart
res@gsnXYBarChartBarWidth 		= 1.0

;res@tmXBOn                		= False                ; turn off tickmarks at bot

res@trYMinF 				= min(var_interp_all)-5.0
res@trYMaxF 				= max(var_interp_all)+5.0
res@trXMinF 				= 0.
res@trXMaxF 				= 50.0

res@tiMainString          		= "Surface heat budget"

x1 					= fspan(2.0,40.0,dimsizes(var2D))
x2 					= fspan(3.0,41.0,dimsizes(var2D))
x3 					= fspan(4.0,42.0,dimsizes(var2D))

res@tmXBMode   				= "Explicit"
res@tmXBValues 				= x2
res@tmXBLabels 				= var2D

res@tmXBLabelAngleF 			= 45
res@tmXBLabelFontHeightF 		= 0.01


res@gsnXYBarChartColors2 		= (/"grey"/)
plot1 					= gsn_csm_xy(wks,x1,var_interp_all(:,0),res)
res@gsnXYBarChartColors2 		= (/"red"/)
plot2 					= gsn_csm_xy(wks,x2,var_interp_all(:,1),res)
res@gsnXYBarChartColors2 		= (/"blue"/)
plot3 					= gsn_csm_xy(wks,x3,var_interp_all(:,2),res)

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
exit


res@trYMinF 				= min(var_interp_diff_all)-5.0 
res@trYMaxF 				= max(var_interp_diff_all)+5.0


res@gsnXYBarChartColors2 		= (/"grey"/)
plot3 					= gsn_csm_xy(wks,x1,var_interp_diff_all(:,0),res)

res@gsnXYBarChartColors2 		= (/"red"/)
plot4 					= gsn_csm_xy(wks,x2,var_interp_diff_all(:,1),res)

res@gsnXYBarChartColors2 		= (/"blue"/)
plot5 					= gsn_csm_xy(wks,x3,var_interp_diff_all(:,2),res)

labels 					= (/"DEF-OBS","NEW-OBS","NEW-DEF"/)
gsn_labelbar_ndc(wks,nboxes,labels,vpx+vph-lbw,0.77,lbres)

draw(plot3)
frame(wks)

exit

if(compare_type(icomp) .eq."amip-obs" .or. compare_type(icomp) .eq."cmip-obs")then
	break
end if


end do; do icomp=
 
                              
end
