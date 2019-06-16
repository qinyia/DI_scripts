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

wks   = gsn_open_wks ("pdf", "./figure/PAI_CMIP_AMIP" ) 

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


result_obs = get_var_2D(fdir(0),fapend(0),casename(0),var(ivar),"obs")
var_obs = result_obs[0]
lat_obs   = result_obs[1]
lon_obs   = result_obs[2]

result_amip = get_var_2D(fdir(1),fapend(1),casename(1),var(ivar),"model")
var_amip = result_amip[0]
lat_amip   = result_amip[1]
lon_amip   = result_amip[2]

result_cmip = get_var_2D(fdir(2),fapend(2),casename(2),var(ivar),"model")
var_cmip = result_cmip[0]
lat_cmip   = result_cmip[1]
lon_cmip   = result_cmip[2]


result = regrid_diff(var_obs,lat_obs,lon_obs,\
	             var_amip,lat_amip,lon_amip,\
		     var_cmip,lat_cmip,lon_cmip)
var_interp = result[0]
var_interp_diff = result[1]

printMinMax(var_interp,False)

; zonal mean
lon1 = 0
lon2 = 360

;; calculate the PAI index
lat_here = var_interp&lat
rad = 4.0*atan(1.0)/180.
clat = cos(lat_here*rad)
clat!0 = "lat"
clat&lat = lat_here

var_interp_NH = wgt_areaave_Wrap(var_interp(:,{0:20},:),clat({0:20}),1.0,1)
var_interp_SH = wgt_areaave_Wrap(var_interp(:,{-20:0},:),clat({-20:0}),1.0,1)
var_interp_tropical = wgt_areaave_Wrap(var_interp(:,{-20:20},:),clat({-20:20}),1.0,1)
var_interp_PAI = (var_interp_NH-var_interp_SH)/var_interp_tropical

;;; this varible is used to put all results from two cases.
if(icase.eq.0)then
	var_interp_NH_all = new((/5/),float)
	var_interp_NH_all(3*icase:3*icase+2) = var_interp_NH

	var_interp_SH_all = new((/5/),float)
	var_interp_SH_all(3*icase:3*icase+2) = var_interp_SH

	var_interp_tropical_all = new((/5/),float)
	var_interp_tropical_all(3*icase:3*icase+2) = var_interp_tropical

	var_interp_PAI_all = new((/5/),float)
	var_interp_PAI_all(3*icase:3*icase+2) = var_interp_PAI
else

var_interp_NH_all(3*icase:3*icase+1) = var_interp_NH(1:2)
var_interp_SH_all(3*icase:3*icase+1) = var_interp_SH(1:2)
var_interp_tropical_all(3*icase:3*icase+1) = var_interp_tropical(1:2)
var_interp_PAI_all(3*icase:3*icase+1) = var_interp_PAI(1:2)

var_interp_PAI_all!0 = "case"
var_interp_PAI_all&case = (/"obs","amip-def","cmip-def","amip-new","cmip-new"/)
end if

delete(var_interp)
delete(var_interp_NH)
delete(var_interp_SH)
delete(var_interp_tropical)
delete(var_interp_PAI)

end do; do icase=

print("cases----------NH-----------SH--------------Tropical----------PAI--------")
print(var_interp_PAI_all&case+"---------"+var_interp_NH_all+"-------"+var_interp_SH_all+"--------"+var_interp_tropical_all+"--------"+var_interp_PAI_all)

AMIP_PAI = new((/3/),float)
AMIP_PAI(0) = var_interp_PAI_all(0)
AMIP_PAI(1) = var_interp_PAI_all(1)
AMIP_PAI(2) = var_interp_PAI_all(3)

CMIP_PAI = new((/3/),float)
CMIP_PAI(0) = var_interp_PAI_all(0)
CMIP_PAI(1) = var_interp_PAI_all(2)
CMIP_PAI(2) = var_interp_PAI_all(4)



; plotting parameters
;*************************
cases = (/"OBS","RH","PDF"/)
colors = (/"black","red","blue"/)
markers = (/16,16,16/)
MarkLineModes = (/"Markers","Markers","Markers"/)
dashes = (/0,2,4/)

;drawNDCGrid(wks)

  res                      = True                 ; plot mods desired

  res@gsnFrame             = False
  res@gsnDraw              = False

  res@vpXF = 0.2
  res@vpYF = 0.85
  res@vpWidthF  = 0.5
  res@vpHeightF = 0.5

  res@xyMarkLineMode = "Markers"               ; Markers *and* lines

  res@tiYAxisString = "PAI_Coupled"
  res@tiXAxisString = "PAI_AMIP"

  res@trYMinF = 0.0
  res@trYMaxF = 0.30
  res@trXMinF = 0.0
  res@trXMaxF = 0.30

  plot(is,0)  = gsn_csm_xy(wks,AMIP_PAI, CMIP_PAI, res )   ; plaace holder

  do i=0,dimsizes(AMIP_PAI)-1

  polyres = True
  polyres@gsMarkerIndex = 16
  polyres@gsMarkerSizeF = 15.0
  polyres@gsMarkerColor = colors(i)

  gsn_polymarker(wks,plot(is,0),AMIP_PAI(i),CMIP_PAI(i),polyres)
 
  end do

  xleg = (/0.20,0.40,0.60/)   ; Location of
  xtxt = (/0.26,0.46,0.66/)   ; legend markers
  yleg = (/0.22,0.22,0.22/)   ; and text
  ytxt = (/0.22,0.22,0.22/)   ; strings.



  gsres = True
  gsres@gsMarkerIndex = 16
  gsres@gsMarkerSizeF = 15.0

  txres = True
  txres@txFontHeightF = 0.020

  do i=0,dimsizes(AMIP_PAI)-1

  gsres@gsMarkerColor = colors(i)

  gsn_polymarker_ndc(wks, xleg(i),yleg(i),gsres)
  gsn_text_ndc(wks, cases(i),xtxt(i),ytxt(i),txres)

  end do

draw(plot)
frame(wks)
exit


;res = get_res_zonal_mean(lat_min,lat_max)

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
