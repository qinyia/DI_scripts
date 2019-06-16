
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

function regrid_diff(var_obs[*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                     var_def[*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		     var_new[*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

;horizontal interpolation
var_interp = lat_lon_interpolate_3vars(var_obs,lat_obs,lon_obs, \
				       var_def,lat_def,lon_def, \
				       var_new,lat_new,lon_new)
var_interp_diff = var_interp

var_interp_diff = var_interp
var_interp_diff(0,:,:) = var_interp(1,:,:) - var_interp(0,:,:) ; def-obs
var_interp_diff(1,:,:) = var_interp(2,:,:) - var_interp(0,:,:) ; new-obs
var_interp_diff(2,:,:) = var_interp(2,:,:) - var_interp(1,:,:) ; new-def

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin


;var=(/"SST","PRECT","CLDTOT"/)
;var=(/"SST","FSNS","FLNS","LHFLX","SHFLX","NetSfcHeat","CLDTOT","PRECT"/)
var=(/"SST","FSNS","FLNS","LHFLX","SHFLX","NetSfcHeat","PRECT"/)

;var=(/"FLNS"/)
;var=(/"LHFLX"/)
;var=(/"CLDTOT"/)

;compare_type = (/"amip-obs","cmip-obs","cmip-amip"/)
;compare_type = (/"cmip-amip"/)
compare_type = (/"cmip-obs"/)
;compare_type = (/"amip-obs"/)

case=(/"def","new"/)

Dir="/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

do icase=0,dimsizes(case)-1
print(case(icase))

do ivar=0,dimsizes(var)-1
print(var(ivar))


	if(compare_type(icomp) .eq."amip-obs")then
		if(var(ivar).eq."SST")then
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
	if(compare_type(icomp) .eq."cmip-amip")then
		if(var(ivar).eq."SST")then
			print("cannot process SST in this compare_type:"+compare_type(icomp))
			exit
		end if
		fdir = (/Dir+"obs/",Dir+"amip/",Dir+"atm/"/)
		if(case(icase).eq."def")then
			casename = (/"","FC5_f09f09_MG15_COSP","BC5_f19g16"/)
		end if
		if(case(icase).eq."new")then
			casename = (/"","new-cldfrac16-1d-MG15-revise-2-COSP","BC5_f19g16_mac2"/)
		end if
	end if



fapend = (/"_ANN_climo.nc", \
           "_ANN_climo.nc", \
	   "_ANN_climo.nc" /)


if(var(ivar).eq."SST")then
	casename(0) = "HadISST_PD"
	fdir = (/Dir+"obs/", \
	         Dir+"ocn/", \
		 Dir+"ocn/"/)
end if

if(var(ivar).eq."PRECT")then
	casename(0) = "GPCP"
end if

if(var(ivar).eq."CLDLOW".or.var(ivar).eq."CLDMED".or.var(ivar).eq."CLDHGH".or.var(ivar).eq."CLDTOT")then
	casename(0) = "CLOUDSAT"
end if

if(var(ivar).eq."FSNS" .or. var(ivar).eq."FLNS" .or. var(ivar).eq."LHFLX" .or. var(ivar).eq."SHFLX".or.var(ivar).eq."NetSfcHeat")then
	casename(0) = "JRA25"
end if



result_obs	= get_var_2D(fdir(0),fapend(0),casename(0),var(ivar),"obs")
var_obs		= result_obs[0]
lat_obs		= result_obs[1]
lon_obs		= result_obs[2]

result_def	= get_var_2D(fdir(1),fapend(1),casename(1),var(ivar),"model")
var_def		= result_def[0]
lat_def		= result_def[1]
lon_def		= result_def[2]

result_new	= get_var_2D(fdir(2),fapend(2),casename(2),var(ivar),"model")
var_new		= result_new[0]
lat_new		= result_new[1]
lon_new		= result_new[2]

printVarSummary(var_obs)

result = regrid_diff(var_obs,lat_obs,lon_obs,\
                     var_def,lat_def,lon_def,\
		     var_new,lat_new,lon_new)

var_interp = result[0]
var_interp_diff = result[1]

;;;;; specific region mean; [0-20S, 180E-140W]

lat1 = -20
lat2 = 0

lon1=180
lon2=220

lat_here = var_interp&lat
rad = 4.0*atan(1.0)/180.
clat = cos(lat_here*rad)
clat!0 = "lat"
clat&lat = lat_here


; area mean over this specific region 
var_interp_areamean = wgt_areaave_Wrap(var_interp(:,{lat1:lat2},{lon1:lon2}),clat({lat1:lat2}),1.0,1) ; obs, def, new
var_interp_diff_areamean = wgt_areaave_Wrap(var_interp_diff(:,{lat1:lat2},{lon1:lon2}),clat({lat1:lat2}),1.0,1) ; def-obs,new-obs,new-def

printVarSummary(var_interp_diff_areamean)

delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/var_interp,var_interp_diff/])
delete(fapend)
delete(fdir)
delete([/lat_here,clat/])

if(ivar.eq.0)then
	var_interp_all = new((/dimsizes(var),3/),float)
	var_interp_all!0 = "vars"
	var_interp_all&vars = var
	var_interp_all!1 = "cases"
	var_interp_all&cases = (/"obs","def","new"/)
	var_interp_all(ivar,:) = (/var_interp_areamean(:)/)

	var_interp_diff_all = new((/dimsizes(var),3/),float)
	var_interp_diff_all!0 = "vars"
	var_interp_diff_all&vars = var
	var_interp_diff_all!1 = "cases"
	var_interp_diff_all&cases = (/"obs","def","new"/)
	var_interp_diff_all(ivar,:) = (/var_interp_diff_areamean(:)/)

else
	var_interp_all(ivar,:) = (/var_interp_areamean(:)/)
	var_interp_diff_all(ivar,:) = (/var_interp_diff_areamean(:)/)
end if


end do ; do ivar=

printVarSummary(var_interp_all)
print(var_interp_all)

;==============================================================================================================
wks = gsn_open_wks("x11","./figure/surface_energy_budget")

cmap1 = read_colormap_file("WhiteBlueGreenYellowRed")
cmap2 = read_colormap_file("ViBlGrWhYeOrRe")
;cmap2 = read_colormap_file("BlueRed")

;plot = new(2,graphic)

res                       = True                 ; plot mods desired         
res@gsnFrame              = False                ; don't advance frame yet

res@tmXMajorGrid 		= True
res@tmXMajorGridThicknessF 	= 1.0
res@tmXMajorGridLineDashPattern	= 2

res@tmYMajorGrid 		= True
res@tmYMajorGridThicknessF 	= 1.0
res@tmYMajorGridLineDashPattern	= 2

res@gsnYRefLine = 0.0


res@gsnXYBarChart         = True                 ; turn on bar chart
res@gsnXYBarChartBarWidth = 1.0

;res@tmXBOn                = False                ; turn off tickmarks at bot

res@trYMinF = min(var_interp_all)-5.0
res@trYMaxF = max(var_interp_all)+5.0
res@trXMinF = 0.
res@trXMaxF = 50.0

res@tiMainString          = "Surface heat budget"

x1 = fspan(2.0,40.0,dimsizes(var))
x2 = fspan(3.0,41.0,dimsizes(var))
x3 = fspan(4.0,42.0,dimsizes(var))

res@tmXBMode = "Explicit"
res@tmXBValues = x2
res@tmXBLabels = var

res@tmXBLabelAngleF = 45
res@tmXBLabelFontHeightF = 0.01


res@gsnXYBarChartColors2 = (/"grey"/)
plot1 = gsn_csm_xy(wks,x1,var_interp_all(:,0),res)
res@gsnXYBarChartColors2 = (/"red"/)
plot2 = gsn_csm_xy(wks,x2,var_interp_all(:,1),res)
res@gsnXYBarChartColors2 = (/"blue"/)
plot3 = gsn_csm_xy(wks,x3,var_interp_all(:,2),res)

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

labels = (/"OBS","DEF","NEW"/)
   
  lbres                    = True          ; labelbar only resources
  lbres@vpWidthF           = 0.2 * vpw     ; labelbar width
  lbres@vpHeightF          = 0.25 * vph     ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.75          ; puts space between color boxes
  lbres@lbFillColors       = (/"grey","red","blue"/) ; labelbar colors
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.008         ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
   
  gsn_labelbar_ndc(wks,nboxes,labels,vpx+vph-lbw,0.77,lbres)
   
draw(plot1) 
frame(wks)


res@trYMinF = min(var_interp_diff_all)-5.0 
res@trYMaxF = max(var_interp_diff_all)+5.0


res@gsnXYBarChartColors2 = (/"grey"/)
plot3 = gsn_csm_xy(wks,x1,var_interp_diff_all(:,0),res)
res@gsnXYBarChartColors2 = (/"red"/)
plot4 = gsn_csm_xy(wks,x2,var_interp_diff_all(:,1),res)
res@gsnXYBarChartColors2 = (/"blue"/)
plot5 = gsn_csm_xy(wks,x3,var_interp_diff_all(:,2),res)

labels = (/"DEF-OBS","NEW-OBS","NEW-DEF"/)
gsn_labelbar_ndc(wks,nboxes,labels,vpx+vph-lbw,0.77,lbres)

draw(plot3)
frame(wks)

exit

if(compare_type(icomp) .eq."amip-obs" .or. compare_type(icomp) .eq."cmip-obs")then
	break
end if

end do; do icase=

end do; do icomp=
 
                              
end
