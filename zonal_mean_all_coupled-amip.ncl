load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff(var_obs[*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                     var_amip[*][*]:numeric,lat_amip[*]:numeric,lon_amip[*]:numeric,\
		     var_cmip[*][*]:numeric,lat_cmip[*]:numeric,lon_cmip[*]:numeric)
begin
;horizontal interpolation
var_interp = lat_lon_interpolate(var_obs,lat_obs,lon_obs, \
				 var_amip,lat_amip,lon_amip, \
				 var_cmip,lat_cmip,lon_cmip)

var_interp_diff = var_interp
var_interp_diff(0,:,:) = var_interp(1,:,:) - var_interp(0,:,:) ; def-obs
var_interp_diff(1,:,:) = var_interp(2,:,:) - var_interp(0,:,:) ; new-obs
var_interp_diff(2,:,:) = var_interp(2,:,:) - var_interp(1,:,:) ; new-obs

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin
;*********************************************
; read in data
;*********************************************

case=(/"def","new"/)

;season=(/"ANN","DJF","MAM","JJA","SON"/)
season=(/"ANN"/)

;var=(/"PRECT","CLDTOT","CLDHGH","CLDMED","CLDLOW"/)
var=(/"PRECT"/)

fdir = (/"/home/share3/lyl/work3/qinyi/double-ITCZ/data/obs/", \
         "/home/share3/lyl/work3/qinyi/double-ITCZ/data/amip/", \
	 "/home/share3/lyl/work3/qinyi/double-ITCZ/data/atm/"/)

do ivar=0,dimsizes(var)-1

wks   = gsn_open_wks ("x11", "./figure/zonal_mean-"+var(ivar)+"-obs-allseason" ) 

plot = new((/dimsizes(season),6/),graphic)

do is = 0,dimsizes(season)-1

fapend = (/"_"+season(is)+"_climo.nc", \
           "_"+season(is)+"_climo.nc", \
	   "_"+season(is)+"_climo.nc" /)

do icase=0,dimsizes(case)-1
print(case(icase))

if(case(icase).eq."def")then
	casename=(/"GPCP","FC5_f19f19_MG15_cosp","BC5_f19g16"/)
end if

if(case(icase).eq."new")then
	casename=(/"GPCP","FC5_f19f19_MG15_mac3_cosp","BC5_f19g16_mac2"/)
end if


if(var(ivar).eq."SST")then
	casename(0) = "HadISST_PD"
end if

if(var(ivar).eq."PRECT")then
	casename(0) = "GPCP"
end if

if(var(ivar).eq."CLDLOW".or.var(ivar).eq."CLDMED".or.var(ivar).eq."CLDHGH".or.var(ivar).eq."CLDTOT")then
	casename(0) = "CLOUDSAT"
end if


if(var(ivar).eq."SST")then
result_obs = get_SST_2D(fdir(0),fapend(0),casename(0),"obs")
var_obs = result_obs[0]
lat_obs   = result_obs[1]
lon_obs   = result_obs[2]

result_amip = get_SST_2D(fdir(1),fapend(1),casename(1),"model")
var_amip = result_amip[0]
lat_amip   = result_amip[1]
lon_amip   = result_amip[2]

result_cmip = get_SST_2D(fdir(2),fapend(2),casename(2),"model")
var_cmip = result_cmip[0]
lat_cmip   = result_cmip[1]
lon_cmip   = result_cmip[2]
end if

if(var(ivar).eq."PRECT")then
result_obs = get_PRECT_2D(fdir(0),fapend(0),casename(0),"obs")
var_obs   = result_obs[0]
lat_obs   = result_obs[1]
lon_obs   = result_obs[2]

result_amip = get_PRECT_2D(fdir(1),fapend(1),casename(1),"model")
var_amip   = result_amip[0]
lat_amip   = result_amip[1]
lon_amip   = result_amip[2]

result_cmip = get_PRECT_2D(fdir(2),fapend(2),casename(2),"model")
var_cmip   = result_cmip[0]
lat_cmip   = result_cmip[1]
lon_cmip   = result_cmip[2]

end if

if(var(ivar).eq."CLDTOT")then
result_obs = get_CLDTOT_2D(fdir(0),fapend(0),casename(0),"obs")
var_obs   = result_obs[0]
lat_obs   = result_obs[1]
lon_obs   = result_obs[2]

result_amip = get_CLDTOT_2D(fdir(1),fapend(1),casename(1),"model")
var_amip   = result_amip[0]
lat_amip   = result_amip[1]
lon_amip   = result_amip[2]

result_cmip = get_CLDTOT_2D(fdir(2),fapend(2),casename(2),"model")
var_cmip   = result_cmip[0]
lat_cmip   = result_cmip[1]
lon_cmip   = result_cmip[2]

end if

if(var(ivar).eq."CLDLOW")then
result_obs = get_CLDLOW_2D(fdir(0),fapend(0),casename(0),"obs")
var_obs   = result_obs[0]
lat_obs   = result_obs[1]
lon_obs   = result_obs[2]

result_amip = get_CLDLOW_2D(fdir(1),fapend(1),casename(1),"model")
var_amip   = result_amip[0]
lat_amip   = result_amip[1]
lon_amip   = result_amip[2]

result_cmip = get_CLDLOW_2D(fdir(2),fapend(2),casename(2),"model")
var_cmip   = result_cmip[0]
lat_cmip   = result_cmip[1]
lon_cmip   = result_cmip[2]

end if

if(var(ivar).eq."CLDMED")then
result_obs = get_CLDMED_2D(fdir(0),fapend(0),casename(0),"obs")
var_obs   = result_obs[0]
lat_obs   = result_obs[1]
lon_obs   = result_obs[2]

result_amip = get_CLDMED_2D(fdir(1),fapend(1),casename(1),"model")
var_amip   = result_amip[0]
lat_amip   = result_amip[1]
lon_amip   = result_amip[2]

result_cmip = get_CLDMED_2D(fdir(2),fapend(2),casename(2),"model")
var_cmip   = result_cmip[0]
lat_cmip   = result_cmip[1]
lon_cmip   = result_cmip[2]

end if

if(var(ivar).eq."CLDHGH")then
result_obs = get_CLDHGH_2D(fdir(0),fapend(0),casename(0),"obs")
var_obs   = result_obs[0]
lat_obs   = result_obs[1]
lon_obs   = result_obs[2]

result_amip = get_CLDHGH_2D(fdir(1),fapend(1),casename(1),"model")
var_amip   = result_amip[0]
lat_amip   = result_amip[1]
lon_amip   = result_amip[2]

result_cmip = get_CLDHGH_2D(fdir(2),fapend(2),casename(2),"model")
var_cmip   = result_cmip[0]
lat_cmip   = result_cmip[1]
lon_cmip   = result_cmip[2]

end if


result = regrid_diff(var_obs,lat_obs,lon_obs,\
	             var_amip,lat_amip,lon_amip,\
		     var_cmip,lat_cmip,lon_cmip)
var_interp = result[0]
var_interp_diff = result[1]

printMinMax(var_interp,False)

; zonal mean
lon1 = 0
lon2 = 360

;printVarSummary(var_interp)

var_interp_zm = dim_avg_n_Wrap(var_interp(:,:,{lon1:lon2}),2)
var_interp_diff_zm = dim_avg_n_Wrap(var_interp_diff(:,:,{lon1:lon2}),2)

var_interp_zm_rev = var_interp_zm(case|:,lat|:)
var_interp_diff_zm_rev = var_interp_diff_zm(case|:,lat|:)

;printVarSummary(var_interp_zm_rev)

;;; this varible is used to put all results from two cases.
if(icase.eq.0)then
var_interp_zm_rev_all = new((/5,dimsizes(var_interp_zm_rev&lat)/),float)
var_interp_zm_rev_all(3*icase:3*icase+2,:) = var_interp_zm_rev
else
var_interp_zm_rev_all(3*icase:3*icase+1,:) = var_interp_zm_rev(1:2,:)
end if

delete(var_interp)
delete(var_interp_zm)
delete(var_interp_zm_rev)

end do; do icase=

; plotting parameters
;*************************

lat_min = -45.
lat_max = 45.

res = get_res_zonal_mean(lat_min,lat_max)

res@gsnCenterString = season(is)
res@gsnRightString = var(ivar)

res@xyMarkLineModes="Lines"
res@xyMarkers=16
res@xyLineThicknessF = 4.0
res@xyDashPatterns = (/0,0,0,1,1/)
res@xyLineColors = (/"grey","red","blue","red","blue"/)

plot(is,0) = gsn_csm_xy(wks,var_interp_zm_rev_all&lat,var_interp_zm_rev_all,res)

;========================
;Attach a legend
lgres = True
lgres@lgLineColors = res@xyLineColors
lgres@lgLineThicknessF = res@xyLineThicknessF
lgres@lgDashIndexes = res@xyDashPatterns
lgres@lgItemTypes = res@xyMarkLineModes
lgres@lgLabelFontHeightF = .10
lgres@vpWidthF = 0.20
lgres@vpHeightF = 0.22
;lgres@lgPerimThicknessF = 1.0
;lgres@lgPerimOn = False

labels = (/"OBS","AMIP_Def","Coupled_Def","AMIP_New","Coupled_New"/)

legend = gsn_create_legend(wks,5,labels,lgres)

amres = True
amres@amJust = "TopLeft"
amres@amParallelPosF = -0.5
amres@amOrthogonalPosF = -0.5 
annoid = gsn_add_annotation(plot(is,0),legend,amres)

draw(plot(is,0))
frame(wks)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;gsn_panel(wks,plot(is,3:5),(/3,1/),False)

;gsn_panel(wks,plot(is,:),(/2,3/),False)

end do ; do is=

;gsn_panel(wks,ndtooned(plot),(/5,6/),False)

delete(var_obs)
delete(lat_obs)
delete(lon_obs)

delete(var_amip)
delete(lat_amip)
delete(lon_amip)

delete(var_cmip)
delete(lat_cmip)
delete(lon_cmip)


delete(plot)
delete(fapend)


end do; do ivar=

                               
end
