load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff(var_1[*][*]:numeric,lat_1[*]:numeric,lon_1[*]:numeric,\
		     var_2[*][*]:numeric,lat_2[*]:numeric,lon_2[*]:numeric)
begin
;horizontal interpolation
var_interp = lat_lon_interpolate_2vars(var_1,lat_1,lon_1, \
				       var_2,lat_2,lon_2)

var_interp_diff = var_interp(0,:,:)
var_interp_diff = var_interp(1,:,:) - var_interp(0,:,:)

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

;var=(/"SST","FSNS","FSNTOA"/)
;var=(/"PRECT","SST","NetSfcHeat","NetSfcRad"/)
var=(/"PRECT","SST","NetSfcHeat","NetSfcRad"/)
;var=(/"SST","NetSfcHeat","NetSfcRad"/)
;var=(/"PRECT","NetSfcHeat"/)

;var=(/"SST"/)

plot = new((/dimsizes(season),6/),graphic)

do is = 0,dimsizes(season)-1

fapend = (/"_"+season(is)+"_climo.nc", \
	   "_"+season(is)+"_climo.nc" /)


do ivar=0,dimsizes(var)-1

wks   = gsn_open_wks ("pdf", "./figure/PAI_change-asymmetry of others" ) 
;wks   = gsn_open_wks ("pdf", "./figure/asymmetry of SST-asymmetry of others" ) 

if(var(ivar).eq."PRECT")then
fdir = (/"/home/share3/lyl/work3/qinyi/double-ITCZ/data/amip/", \
	 "/home/share3/lyl/work3/qinyi/double-ITCZ/data/atm/"/)
end if

if(var(ivar).eq."SST")then
fdir = (/"/home/share3/lyl/work3/qinyi/double-ITCZ/data/obs/", \
	 "/home/share3/lyl/work3/qinyi/double-ITCZ/data/ocn/"/)
end if

if(var(ivar).ne."SST" .and. var(ivar).ne."PRECT")then
fdir = (/"/home/share3/lyl/work3/qinyi/double-ITCZ/data/obs/", \
	 "/home/share3/lyl/work3/qinyi/double-ITCZ/data/amip/"/)
end if


do icase=0,dimsizes(case)-1

; note: for total precipitation rate, it is the difference of PAI between Coupled and AMIP --> get PAI index change from AMIP to CMIP.
; for surface temperature bias, it is the difference of SST bias (Coupled-OBS)  between NH and SH --> asymmetry of surf temp bias.
; for other variables, it is the difference of "variable" in AMIP between NH and SH, not its bias. --> asymmetry of "variable".
; the basic logic here is: asymmetry of "other variable" in AMIP--> asymmetry of SST bias in coupled simulation --> PAI change from AMIP to CMIP.

if(case(icase).eq."def")then
	casename=(/"","FC5_f19f19_MG15_cosp"/)
end if

if(case(icase).eq."new")then
	casename=(/"","FC5_f19f19_MG15_mac3_cosp"/)
end if

if(var(ivar).eq."PRECT")then ; coupled and AMIP
	if(case(icase).eq."def")then
		casename=(/"FC5_f19f19_MG15_cosp","BC5_f19g16"/)
	end if
	
	if(case(icase).eq."new")then
		casename=(/"FC5_f19f19_MG15_mac3_cosp","BC5_f19g16_mac2"/)
	end if
end if


if(var(ivar).eq."SST")then ; coupled and OBS (AMIP)
	if(case(icase).eq."def")then
		casename=(/"HadISST_PD","BC5_f19g16"/)
	end if
	
	if(case(icase).eq."new")then
		casename=(/"HadISST_PD","BC5_f19g16_mac2"/)
	end if
end if

; other variables: only AMIP
if(var(ivar).eq."FSNS".or.var(ivar).eq."FLNS".or.\
   var(ivar).eq."LHFLX".or.var(ivar).eq."SHFLX".or.\
   var(ivar).eq."NetSfcHeat".or.var(ivar).eq."NetSfcRad")then

	casename(0) = "JRA25"
end if

if(var(ivar).eq."FSNTOA")then
	casename(0) = "CERES-EBAF"
end if

if(var(ivar).eq."PRECT")then
result_obs = get_var_2D(fdir(0),fapend(0),casename(0),var(ivar),"model")
var_obs   = result_obs[0]
lat_obs   = result_obs[1]
lon_obs   = result_obs[2]

result_cmip = get_var_2D(fdir(1),fapend(1),casename(1),var(ivar),"model")
var_cmip   = result_cmip[0]
lat_cmip   = result_cmip[1]
lon_cmip   = result_cmip[2]

else

result_obs = get_var_2D(fdir(0),fapend(0),casename(0),var(ivar),"obs")
var_obs   = result_obs[0]
lat_obs   = result_obs[1]
lon_obs   = result_obs[2]

result_cmip = get_var_2D(fdir(1),fapend(1),casename(1),var(ivar),"model")
var_cmip   = result_cmip[0]
lat_cmip   = result_cmip[1]
lon_cmip   = result_cmip[2]


end if


result = regrid_diff(var_obs,lat_obs,lon_obs,\
		     var_cmip,lat_cmip,lon_cmip)
var_interp = result[0]
var_interp_diff = result[1]

; zonal mean
lat1 = -20
lat2 = 0
lat3 = 20


;; calculate the asy index
lat_here = var_interp&lat
rad = 4.0*atan(1.0)/180.
clat = cos(lat_here*rad)
clat!0 = "lat"
clat&lat = lat_here

if(var(ivar).eq."PRECT")then ; for PAI, get the PAI for coupled and AMIP, respectively.
	var_interp_NH = wgt_areaave_Wrap(var_interp(:,{lat2:lat3},:),clat({lat2:lat3}),1.0,1)
	var_interp_SH = wgt_areaave_Wrap(var_interp(:,{lat1:lat2},:),clat({lat1:lat2}),1.0,1)
	var_interp_tropical = wgt_areaave_Wrap(var_interp(:,{lat1:lat3},:),clat({lat1:lat3}),1.0,1)
	var_interp_asy = (var_interp_NH-var_interp_SH)/var_interp_tropical ; PAI index
end if

if(var(ivar).eq."SST")then ; asymmetry of surf temp bias. use the difference b/t coupled and obs to get the asymmetry directly.
	var_interp_NH = wgt_areaave_Wrap(var_interp_diff({lat2:lat3},:),clat({lat2:lat3}),1.0,1)
	var_interp_SH = wgt_areaave_Wrap(var_interp_diff({lat1:lat2},:),clat({lat1:lat2}),1.0,1)
	var_interp_asy = var_interp_NH-var_interp_SH ; asymmetry of SST bias (CMIP-OBS)
end if

if(var(ivar).ne."PRECT" .and. var(ivar).ne."SST")then ; asymmetry of net sfc heat flux/TOA SW radiation, et al. Actually, later we just need to AMIP result. but obs result is used as a comparison.
	var_interp_NH = wgt_areaave_Wrap(var_interp(:,{lat2:lat3},:),clat({lat2:lat3}),1.0,1)
	var_interp_SH = wgt_areaave_Wrap(var_interp(:,{lat1:lat2},:),clat({lat1:lat2}),1.0,1)
	var_interp_asy = var_interp_NH-var_interp_SH ; asymmetry of "other vars" in obs. and AMIP simulation
end if


;;; "xxx_all": this varible is used to put all results from two cases and obs.
if(icase.eq.0)then
	var_interp_NH_all = new((/3/),float) ; obs, def, new
	if(var(ivar).eq."SST")then 
		var_interp_NH_all(icase) = 0.	; in obs, no surf temp bias for SST.
		var_interp_NH_all(icase+1) = var_interp_NH ; in new, SST bias b/t coupled and obs in NH.
	end if
	if(var(ivar).eq."PRECT")then
		var_interp_NH_all(icase) = 0. ; in obs, not caculate the pr. change because there is only obs and no comparison.
		var_interp_NH_all(icase+1) =(/var_interp_NH(1)-var_interp_NH(0)/) ; in NH, pr. change from coupled to AMIP in def simulation.
	end if

	if(var(ivar).ne."SST" .and. var(ivar).ne."PRECT")then
		var_interp_NH_all(icase:icase+1) = var_interp_NH(:) ; in NH, "vars" in obs. and def simulation from AMIP run.
	end if

	var_interp_SH_all = new((/3/),float)
	if(var(ivar).eq."SST")then
		var_interp_SH_all(icase) = 0.
		var_interp_SH_all(icase+1) = var_interp_SH ; SST bias in SH
	end if
	if(var(ivar).eq."PRECT")then
		var_interp_SH_all(icase) = 0.
	        var_interp_SH_all(icase+1) =(/var_interp_SH(1)-var_interp_SH(0)/) ; pr .change from coupled to AMIP in SH
	end if
	if(var(ivar).ne."SST" .and. var(ivar).ne."PRECT")then
		var_interp_SH_all(icase:icase+1) = var_interp_SH(:) ; in SH, "vars" in obs and AMIP run.
	end if


	var_interp_asy_all = new((/3/),float)
	if(var(ivar).eq."SST")then
		var_interp_asy_all(icase) = 0. ; no asymmetry of SST bias in obs
		var_interp_asy_all(icase+1) = var_interp_asy ; asymmetry of SST bias in different cases.
	end if
	if(var(ivar).eq."PRECT")then
		var_interp_asy_all(icase) = 0. ; no calculation of PAI difference from coupled to AMIP in OBS
	        var_interp_asy_all(icase+1) =(/var_interp_asy(1)-var_interp_asy(0)/) ; change of PAI from coupled to AMIP
	end if
	if(var(ivar).ne."SST" .and. var(ivar).ne."PRECT")then
		var_interp_asy_all(icase:icase+1) = (/var_interp_asy(:)/) ; asymmetry (NH-SH) of "other vars" in AMIP
	end if



else
	if(var(ivar).eq."SST")then
	      var_interp_NH_all(icase+1) = var_interp_NH
	      var_interp_SH_all(icase+1) = var_interp_SH
	     var_interp_asy_all(icase+1) = var_interp_asy
	end if
	if(var(ivar).eq."PRECT")then
	      var_interp_NH_all(icase+1) =(/var_interp_NH(1)-var_interp_NH(0)/) ; new
	      var_interp_SH_all(icase+1) =(/var_interp_SH(1)-var_interp_SH(0)/)
	     var_interp_asy_all(icase+1) =(/var_interp_asy(1)-var_interp_asy(0)/)
	end if

	if(var(ivar).ne."SST" .and. var(ivar).ne."PRECT")then
	      var_interp_NH_all(icase+1) =(/var_interp_NH(1)/) ; new
	      var_interp_SH_all(icase+1) =(/var_interp_SH(1)/)
	     var_interp_asy_all(icase+1) =(/var_interp_asy(1)/)
	end if


var_interp_asy_all!0 = "case"
var_interp_asy_all&case = (/"OBS","RH","PDF"/)
end if


delete(lat_here)
delete(clat)
delete(var_interp)
delete(var_interp_diff)
delete(var_interp_NH)
delete(var_interp_SH)
delete(var_interp_asy)

delete(var_obs)
delete(lat_obs)
delete(lon_obs)

delete(var_cmip)
delete(lat_cmip)
delete(lon_cmip)

end do; do icase=


print("cases----------NH-----------SH----------asy--------")
print(var_interp_asy_all&case+"		"+var_interp_NH_all+"		"+var_interp_SH_all+"		"+var_interp_asy_all)

;;; output to one table file.
header = (/"--------------------------------------------",\
	   var(ivar),\
           "cases--------NH--------SH--------asy--------",\
	   "--------------------------------------------"/)

hlist = [/header/]

alist = [/var_interp_asy_all&case, var_interp_NH_all, var_interp_SH_all, var_interp_asy_all/]

fname_t = "table"+var(ivar)+".txt"
write_table(fname_t,"w",hlist, "%s")
write_table(fname_t,"a",alist, "%s%10.3f%10.3f%10.3f")

if(ivar.eq.0)then
asy_multi_vars = new((/dimsizes(var),dimsizes(var_interp_asy_all)/),float) ; number of vars, cases
end if

asy_multi_vars(ivar,:) = var_interp_asy_all

print(asy_multi_vars)
end do; do ivar=

; plotting parameters
;*************************
cases = (/"OBS","RH","PDF"/)
colors = (/"black","red","blue"/)
markers = (/16,16,16/)
MarkLineModes = (/"Markers","Markers","Markers"/)
dashes = (/0,2,4/)


res                      = True                 ; plot mods desired

res@gsnFrame             = False
res@gsnDraw              = False

res@vpXF = 0.2
res@vpYF = 0.85
res@vpWidthF  = 0.5
res@vpHeightF = 0.5

res@xyMarkLineMode = "Markers"               ; Markers *and* lines

do ivar2 = 1,dimsizes(var)-1

;drawNDCGrid(wks)


if(var(0).eq."SST")then
	res@tiXAxisString = "asymmetry of "+var(ivar2)+" in AMIP"
	res@tiYAxisString = "asymmetry of "+var(0)+" bias"
end if
if(var(0).eq."PRECT")then
	res@tiYAxisString = "PAI change (Coupled-AMIP)"
	if(var(ivar2).eq."SST")then
		res@tiXAxisString = "asymmetry of "+var(ivar2)+" bias"
	else
		res@tiXAxisString = "asymmetry of "+var(ivar2)+" in AMIP"
	end if
end if


  res@trYMinF = min(asy_multi_vars(0,:))-0.1
  res@trYMaxF = max(asy_multi_vars(0,:))+0.1
  res@trXMinF = min(asy_multi_vars(ivar2,:))-0.5
  res@trXMaxF = max(asy_multi_vars(ivar2,:))+0.5

plot(is,0)  = gsn_csm_xy(wks,asy_multi_vars(ivar2,:), asy_multi_vars(0,:), res )   ; plaace holder


do i=0,dimsizes(asy_multi_vars(0,:))-1

polyres = True
polyres@gsMarkerIndex = 16
polyres@gsMarkerSizeF = 15.0
polyres@gsMarkerColor = colors(i)

gsn_polymarker(wks,plot(is,0),asy_multi_vars(ivar2,i),asy_multi_vars(0,i),polyres)

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

do i=0,dimsizes(asy_multi_vars(0,:))-1

gsres@gsMarkerColor = colors(i)

gsn_polymarker_ndc(wks, xleg(i),yleg(i),gsres)
gsn_text_ndc(wks, cases(i),xtxt(i),ytxt(i),txres)

end do

draw(plot)
frame(wks)

end do ; do ivar2=


;res = get_res_zonal_mean(lat_min,lat_max)


;gsn_panel(wks,ndtooned(plot),(/5,6/),False)

delete(plot)
delete(fapend)


end do ; do is=

                               
end
