load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

; 2017-12-19: it seems that horizontal advection should be vertically integrated, rather than vertically averaged. --> no!
; 2017-12-19: the vertical integration should be problematic. it should be weighted...

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff_TimeDepthLatLon_2var(var_1[*][*][*][*]:numeric,lat_1[*]:numeric,lon_1[*]:numeric,\
		                     var_2[*][*][*][*]:numeric,lat_2[*]:numeric,lon_2[*]:numeric)
begin
var_1!1		= "z_t"

ntime		= dimsizes(var_1&time)
z_t		= var_1&z_t
ndepth		= dimsizes(z_t)
nlat		= dimsizes(var_1&lat)
nlon		= dimsizes(var_1&lon)

var_interp	= new((/2,ntime,ndepth,nlat,nlon/),"float")

var_interp(0,:,:,:,:)		= var_1(:,:,:,:)
var_interp(1,:,:,:,:)		= var_2(:,:,:,:)

var_interp!0 = "case"
var_interp&case = (/0,1/)
var_interp!1 = "time"
var_interp&time = ispan(0,ntime-1,1)
var_interp!2 = "depth"
var_interp&depth = z_t
var_interp!3 = "lat"
var_interp&lat = var_1&lat
var_interp!4 = "lon"
var_interp&lon = var_1&lon


var_interp_diff = var_interp(0,:,:,:,:)
var_interp_diff = var_interp(1,:,:,:,:) - var_interp(0,:,:,:,:)

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff_TimeLatLon_2var(var_1[*][*][*]:numeric,lat_1[*]:numeric,lon_1[*]:numeric,\
		                     var_2[*][*][*]:numeric,lat_2[*]:numeric,lon_2[*]:numeric)
begin

ntime		= dimsizes(var_1&time)
nlat		= dimsizes(var_1&lat)
nlon		= dimsizes(var_1&lon)

var_interp	= new((/2,ntime,nlat,nlon/),"float")

var_interp(0,:,:,:)		= var_1
var_interp(1,:,:,:)		= var_2

var_interp!0 = "case"
var_interp&case = (/0,1/)
var_interp!1 = "time"
var_interp&time = ispan(0,ntime-1,1)
var_interp!2 = "lat"
var_interp&lat = var_1&lat
var_interp!3 = "lon"
var_interp&lon = var_1&lon


var_interp_diff = var_interp(0,:,:,:)
var_interp_diff = var_interp(1,:,:,:) - var_interp(0,:,:,:)

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin

system("date")

;*********************************************
; read in data
;*********************************************

;season=(/"ANN","DJF","MAM","JJA","SON"/)

season			=(/"ANN"/)
varTimeDepthLatLon	=(/"TEMP","UVEL","VVEL","WVEL"/)

varTimeLatLon		=(/"HMXL","SHF"/)

compare_type 		= (/"cmip-obs"/)

case			= (/"def","new"/)

Output_TimeDepthLatLon		= True
Output_TimeLatLon		= True



;Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"
Dir			= "/home/lyl/WORK4/cesm1_2_1/archive/"

compset			= (/"ocn"/)

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

do icompset=0,dimsizes(compset)-1

if(compare_type(icomp) .eq."cmip-obs")then
Dir			= "/home/lyl/WORK1/lyl_cesm1_2_1/archive/"

	casename 	= (/"BC5_f19g16","BC5_f19g16_mac2"/)
	if(compset(icompset).eq."ocn")then
		fdir_append 	= (/"/ocn/hist/","/ocn/hist/"/)
	end if
	if(compset(icompset).eq."atm")then
		fdir_append	= (/"/atm/hist/","/atm/hist/"/)
	end if
end if


;====================================
; deal with TimeDepthLatLon variables
;====================================

do is = 0,dimsizes(season)-1

do ivar=0,dimsizes(varTimeDepthLatLon)-1


fappend 			= (/"_"+season(is)+"_climo.nc", \
          	            	    "_"+season(is)+"_climo.nc", \
	  	            	    "_"+season(is)+"_climo.nc" /) 
;;; for normal variables

result_def		= get_var_ocn_TimeDepthLatLon_addfiles(Dir,fdir_append(0),casename(0),varTimeDepthLatLon(ivar),"model")
var_def	 	       	= result_def[0]
lat_def	 	       	= result_def[1]
lon_def	 	       	= result_def[2]
time_def		= result_def[3]
depth_def		= result_def[4]

result_new		= get_var_ocn_TimeDepthLatLon_addfiles(Dir,fdir_append(1),casename(1),varTimeDepthLatLon(ivar),"model")
var_new			= result_new[0]
lat_new			= result_new[1]
lon_new			= result_new[2]
time_new		= result_new[3]
depth_new		= result_new[4]

print("read TimeDepthLatLon success")
system("date")

result 			= regrid_diff_TimeDepthLatLon_2var(var_def,lat_def,lon_def,\
		             		              var_new,lat_new,lon_new)
var_interp 		= result[0]
var_interp_diff 	= result[1]

delete([/var_def,depth_def,lat_def,lon_def,var_new,lat_new,lon_new,depth_new/])

print("regrid TimeDepthLatLon success")
system("date")

printVarSummary(var_interp)

index			= dimsizes(var_interp)

if(ivar.eq.0)then
	vars_TimeDepthLatLon_tmp 	= new((/dimsizes(varTimeDepthLatLon),index(0),index(1),index(2),index(3),index(4)/),typeof(var_interp))
	vars_TimeDepthLatLon_tmp!0 	= "vars"
	vars_TimeDepthLatLon_tmp&vars 	= varTimeDepthLatLon
end if
vars_TimeDepthLatLon_tmp(ivar,:,:,:,:,:) 	= var_interp

delete(index)
delete(var_interp)
delete(var_interp_diff)

;;========================================================
;; write into one NC file to check them
;;========================================================
;system("/bin/rm -f vars_TimeDepthLatLon_"+varTimeDepthLatLon(ivar)+"_"+season(is)+"_"+compset(icompset)+"_"+compare_type(icomp)+".nc")
;
;setfileoption("nc","Format","LargeFile")
;ncdf = addfile("vars_TimeDepthLatLon_"+varTimeDepthLatLon(ivar)+"_"+season(is)+"_"+compset(icompset)+"_"+compare_type(icomp)+".nc" ,"c")
;
;	ncdf->$varTimeDepthLatLon(ivar)$		= vars_TimeDepthLatLon(ivar,:,:,:,:,:)
;
;print("write nc successfully")


end do; do ivar=


;====================================
; deal with TimeLatLon variables
;====================================

do ivar=0,dimsizes(varTimeLatLon)-1


fappend 			= (/"_"+season(is)+"_climo.nc", \
          	            	    "_"+season(is)+"_climo.nc", \
	  	            	    "_"+season(is)+"_climo.nc" /) 
;;; for normal variables

result_def		= get_var_ocn_TimeLatLon_addfiles(Dir,fdir_append(0),casename(0),varTimeLatLon(ivar),"model")
var_def	 	       	= result_def[0]
lat_def	 	       	= result_def[1]
lon_def	 	       	= result_def[2]
time_def		= result_def[3]

result_new		= get_var_ocn_TimeLatLon_addfiles(Dir,fdir_append(1),casename(1),varTimeLatLon(ivar),"model")
var_new			= result_new[0]
lat_new			= result_new[1]
lon_new			= result_new[2]
time_new		= result_new[3]

print("read TimeLatLon success")
system("date")

result 			= regrid_diff_TimeLatLon_2var(var_def,lat_def,lon_def,\
		             		              var_new,lat_new,lon_new)
var_interp 		= result[0]
var_interp_diff 	= result[1]

delete([/var_def,lat_def,lon_def,var_new,lat_new,lon_new/])

print("regrid TimeLatLon success")
system("date")

printVarSummary(var_interp)

index			= dimsizes(var_interp)

if(ivar.eq.0)then
	vars_TimeLatLon_tmp 	= new((/dimsizes(varTimeLatLon),index(0),index(1),index(2),index(3)/),typeof(var_interp))
	vars_TimeLatLon_tmp!0 	= "vars"
	vars_TimeLatLon_tmp&vars = varTimeLatLon
end if
vars_TimeLatLon_tmp(ivar,:,:,:,:) 	= var_interp

delete(index)
delete(var_interp)
delete(var_interp_diff)

;;========================================================
;; write into one NC file to check them
;;========================================================
;system("/bin/rm -f vars_TimeLatLon_"+varTimeLatLon(ivar)+"_"+season(is)+"_"+compset(icompset)+"_"+compare_type(icomp)+".nc")
;
;setfileoption("nc","Format","LargeFile")
;ncdf = addfile("vars_TimeLatLon_"+varTimeLatLon(ivar)+"_"+season(is)+"_"+compset(icompset)+"_"+compare_type(icomp)+".nc" ,"c")
;
;	ncdf->$varTimeLatLon(ivar)$		= vars_TimeLatLon(ivar,:,:,:,:)
;
;print("write nc successfully")


end do; do ivar=

print("calculation starts here...")
system("date")

;==============================================
; get the depth index of mixed layer depth
; get the averaged values over the mixed layer.
; get the entrainment velocity.
;==============================================

latS			= -20.
latE			= 20.

lonS			= 160.
lonE			= 280.

vars_TimeDepthLatLon	= vars_TimeDepthLatLon_tmp(:,:,:,:,{latS:latE},{lonS:lonE})
vars_TimeLatLon		= vars_TimeLatLon_tmp(:,:,:,{latS:latE},{lonS:lonE})

; (vars,case,time,lat,lon)
HMXL 			= 	vars_TimeLatLon(0,:,:,:,:) ;unit is m.
depth 			= 	vars_TimeDepthLatLon&depth ; unit is cm.
depth			= 	depth/100.

printVarSummary(HMXL)

; (vars,case,time,depth,lat,lon)
TEMP			= 	vars_TimeDepthLatLon(0,:,:,:,:,:)
UVEL			= 	vars_TimeDepthLatLon(1,:,:,:,:,:)
VVEL			= 	vars_TimeDepthLatLon(2,:,:,:,:,:)
WVEL			= 	vars_TimeDepthLatLon(3,:,:,:,:,:)


nlat			= 	dimsizes(vars_TimeLatLon&lat)
nlon			= 	dimsizes(vars_TimeLatLon&lon)
ncase			=	dimsizes(vars_TimeLatLon&case)
ntime			= 	dimsizes(vars_TimeLatLon&time)
ndepth			= 	dimsizes(vars_TimeDepthLatLon&depth)

lat			= 	vars_TimeLatLon&lat
lon			=	vars_TimeLatLon&lon
case			=	vars_TimeLatLon&case
time			= 	vars_TimeLatLon@time


; for T, U, V, need to get the averaged values from sea surface to mixed layer depth.
; for W, should get the entrainment velocity, which is the vertical velocity at the mixed layer base.

index_mld 			= 	new((/ncase,ntime,nlat,nlon/),integer)
index_mld@longname		= 	"index of mixed layer depth in Depth dimension"
index_mld_tmp 			= 	index_mld

vars_TimeDepthLatLon_mldmean	=	new((/dimsizes(varTimeDepthLatLon),ncase,ntime,nlat,nlon/),float)

; for T, U, V, need to get the averaged values from sea surface to mixed layer depth.
; for W, should get the entrainment velocity, which is the vertical velocity at the mixed layer base.
w_mld				=	new((/ncase,ntime,nlat,nlon/),float)

; the temperature difference between the mixed layer and just below the mixed layer
; qinyi 2017-12-12 12:30:18 corrected
; the temperature difference should be the mixed layer mean temperature and the sublayer temperature
deltaT_mld_old			= 	new((/ncase,ntime,nlat,nlon/),float)
deltaT_mld			= 	new((/ncase,ntime,nlat,nlon/),float)

; 2017-12-18 23:07:24 revised
; the calcualtion of u and v advection terms are based on each level and then get the vertical mean.
; so move the calculation of horizontal advection in front of the do loops.
u_dTdx				= 	new((/ncase,ntime,nlat,nlon/),float)
v_dTdy				= 	new((/ncase,ntime,nlat,nlon/),float)

;========================================================
; calculate horizontal advection: zonal and meridional
;========================================================
; get latitudinal and longitudinal gradient
; here, leave one small problem: how to directly use "TEMP" as the index?
; temperature is the first variable in varTimeDepthLatLon, so its index is 0.
; grad_latlon_cfd only can deal with 2-4 dimensions variables. 
; so I have to do loop to calculate two cases' results.

u_dTdx_lev	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
v_dTdy_lev	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
T_grad_dx	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
T_grad_dy	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)

u_dTdx_lev!0			= "case"
u_dTdx_lev&case			= vars_TimeLatLon&case

v_dTdy_lev!0			= "case"
v_dTdy_lev&case			= vars_TimeLatLon&case

T_grad_dx!0			= "case"
T_grad_dx&case			= vars_TimeLatLon&case

T_grad_dy!0			= "case"
T_grad_dy&case			= vars_TimeLatLon&case


; get gradients
do icase=0,ncase-1

; 2017-12-17: problem: T_grad_dx gets some strange values over north pole. needs to be solved.
; 2017-12-18: the problem comes from the regrided data after interpolation: fv1.9x2.5.
; I am still not sure why it has this problem. Now I change the grid into 1x1d before, and that problem disappears.
gradLatLon  			= 	grad_latlon_cfd (TEMP(icase,:,:,:,:), vars_TimeDepthLatLon&lat, vars_TimeDepthLatLon&lon, False, False)
T_grad_dy(icase,:,:,:,:)  	= 	gradLatLon[0] 
T_grad_dx(icase,:,:,:,:)  	= 	gradLatLon[1]

delete(gradLatLon)
; zonal advection: u_dTdx
u_dTdx_lev(icase,:,:,:,:)	= 	T_grad_dx(icase,:,:,:,:)
u_dTdx_lev(icase,:,:,:,:)	= 	UVEL(icase,:,:,:,:)*T_grad_dx(icase,:,:,:,:)

; meridional advection: v_dTdy
v_dTdy_lev(icase,:,:,:,:)	= 	T_grad_dy(icase,:,:,:,:)
v_dTdy_lev(icase,:,:,:,:)	= 	VVEL(icase,:,:,:,:)*T_grad_dy(icase,:,:,:,:)

end do ; do icase=

T_grad_dx@long_name 		= 	"zonal gradient"
T_grad_dy@long_name 		= 	"meridional gradient"
T_grad_dx@units     		= 	"degC/m"
T_grad_dy@units     		= 	"degC/m"

printVarSummary(u_dTdx_lev)
printVarSummary(T_grad_dx)


; Get weights for depths
wgt				= new((/ndepth/),float)
do idepth=0,ndepth-2
	wgt(idepth)		= depth(idepth+1)-depth(idepth)
end do
wgt(ndepth-1)			= depth(ndepth-1)-depth(ndepth-2)

print(wgt)

; calculations on each grid point

do icase= 0,ncase-1
do itime = 0,ntime-1
;do itime = 0,1
do ilat = 0,nlat-1
do ilon = 0,nlon-1

print(icase+"  "+itime+"  "+ilat+"  "+ilon)

index_mld_tmp(icase,itime,ilat,ilon) 			= 	closest_val(HMXL(icase,itime,ilat,ilon),depth)	

if(.not.ismissing(index_mld_tmp(icase,itime,ilat,ilon)).and.index_mld_tmp(icase,itime,ilat,ilon).gt.1)then
	if(depth(index_mld_tmp(icase,itime,ilat,ilon)).gt.HMXL(icase,itime,ilat,ilon))then
		index_mld(icase,itime,ilat,ilon)	=	index_mld_tmp(icase,itime,ilat,ilon)-1
	else
		index_mld(icase,itime,ilat,ilon)	=	index_mld_tmp(icase,itime,ilat,ilon)
	end if
end if

if(.not.ismissing(index_mld(icase,itime,ilat,ilon)))then
	vars_TimeDepthLatLon_mldmean(:,icase,itime,ilat,ilon) 	= dim_avg_wgt_n_Wrap(vars_TimeDepthLatLon(:,icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,1)

	w_mld(icase,itime,ilat,ilon)				= (/vars_TimeDepthLatLon(3,icase,itime,index_mld(icase,itime,ilat,ilon),ilat,ilon)/)

	deltaT_mld(icase,itime,ilat,ilon)			= (/vars_TimeDepthLatLon_mldmean(0,icase,itime,ilat,ilon)-vars_TimeDepthLatLon(0,icase,itime,index_mld(icase,itime,ilat,ilon)+2,ilat,ilon)/)

	; vertical integral of horizontal advection
	; 2017-12-19 revised
	u_dTdx(icase,itime,ilat,ilon)				= dim_avg_wgt_n_Wrap(u_dTdx_lev(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	v_dTdy(icase,itime,ilat,ilon)				= dim_avg_wgt_n_Wrap(v_dTdy_lev(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	
else
	vars_TimeDepthLatLon_mldmean(:,icase,itime,ilat,ilon) 	=	vars_TimeDepthLatLon@_FillValue
	w_mld(icase,itime,ilat,ilon)				= 	vars_TimeDepthLatLon@_FillValue
	deltaT_mld_old(icase,itime,ilat,ilon)			=	vars_TimeDepthLatLon@_FillValue
	deltaT_mld(icase,itime,ilat,ilon)			=	vars_TimeDepthLatLon@_FillValue
	u_dTdx(icase,itime,ilat,ilon)				=	vars_TimeDepthLatLon@_FillValue
	v_dTdy(icase,itime,ilat,ilon)				=	vars_TimeDepthLatLon@_FillValue
end if

end do ; do ilon=
end do ; do ilat=
end do ; do itime=
end do ; do icase=

system("date")


vars_TimeDepthLatLon_mldmean!0			= "vars"
vars_TimeDepthLatLon_mldmean&vars			= vars_TimeDepthLatLon&vars
vars_TimeDepthLatLon_mldmean!1			= "case"
vars_TimeDepthLatLon_mldmean&case			= vars_TimeDepthLatLon&case
vars_TimeDepthLatLon_mldmean!2			= "time"
vars_TimeDepthLatLon_mldmean&time			= vars_TimeDepthLatLon&time
vars_TimeDepthLatLon_mldmean!3			= "lat"
vars_TimeDepthLatLon_mldmean&lat			= vars_TimeDepthLatLon&lat
vars_TimeDepthLatLon_mldmean!4			= "lon"
vars_TimeDepthLatLon_mldmean&lon			= vars_TimeDepthLatLon&lon

copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),w_mld)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),deltaT_mld_old)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),deltaT_mld)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),index_mld)

copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),u_dTdx)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),v_dTdy)

printVarSummary(vars_TimeDepthLatLon_mldmean)
printVarSummary(u_dTdx)
printMinMax(u_dTdx,False)

Temp_mldmean 		= 	vars_TimeDepthLatLon_mldmean(0,:,:,:,:)
u_mldmean		= 	vars_TimeDepthLatLon_mldmean(1,:,:,:,:)
v_mldmean		= 	vars_TimeDepthLatLon_mldmean(2,:,:,:,:)

;========================================================
; calculate entrainment term
;========================================================

entrain			= HMXL
entrain 		= w_mld*deltaT_mld/HMXL
entrain			= entrain*(-1.0)

entrain 		= entrain*30*24*3600.	;change from K/s to K/month

entrain@longname	= "entrainment term"	
entrain@units		= "degC/month"

printVarSummary(entrain)

;========================================================
; calculate horizontal advection: zonal and meridional
;========================================================
u_dTdx			=	u_dTdx*(-1.)
;change unit from K/s to K/month
u_dTdx 			= 	u_dTdx*30*24*3600

u_dTdx@long_name	= 	"advection term: u_dTdx"
u_dTdx@units 		= 	"degC/month"

printVarSummary(u_dTdx)

; meridional advection: v_dTdy
v_dTdy			= 	v_dTdy*(-1.)
;change units
v_dTdy 			= 	v_dTdy*30*24*3600

v_dTdy@longname 	= 	"advection term: v_dTdy"
v_dTdy@units 		= 	"degC/month"

printVarSummary(v_dTdy)

;========================================================
; calculate the net energy flux term
;========================================================

rho_sw 			= 1030.		; kg/m3
cp_sw			= 4.0e3		; J/kg/K

SHF			= vars_TimeLatLon(1,:,:,:,:) ; total surface heat flux: (case,lat,lon)
printVarSummary(SHF)

netFlux 		= SHF/rho_sw/cp_sw/HMXL
copy_VarCoords(SHF,netFlux)

netFlux			= netFlux*30*24*3600.
netFlux@longname 	= "net energy flux term"
netFlux@units		= "degC/month"

printVarSummary(netFlux)

;========================================================
; calculate the sum of four terms:
; net energy flux: netFlux
; zonal advection: u_dTdx
; meridional advection: v_dTdy
; entrainment: entrain
;========================================================

sumall 			= netFlux + u_dTdx + v_dTdy + entrain
copy_VarCoords(netFlux,sumall)
residual		= (-1.)*sumall
copy_VarCoords(netFlux,residual)

sumall@longname 	= "sum of four terms"
sumall@units 		= "degC/month"

residual@longname 	= "residual"
residual@units 		= "degC/month"

printVarSummary(sumall)

;========================================================
; change units from degC/month to W/m2
;========================================================

factor			= 30.*24.*3600.
residual_Wm2		= residual*cp_sw*rho_sw*HMXL/factor
sumall_Wm2		= sumall*cp_sw*rho_sw*HMXL/factor
netFlux_Wm2		= netFlux*cp_sw*rho_sw*HMXL/factor
u_dTdx_Wm2		= u_dTdx*cp_sw*rho_sw*HMXL/factor
v_dTdy_Wm2		= v_dTdy*cp_sw*rho_sw*HMXL/factor
entrain_Wm2		= entrain*cp_sw*rho_sw*HMXL/factor

copy_VarCoords(residual,residual_Wm2)
copy_VarCoords(sumall,sumall_Wm2)
copy_VarCoords(netFlux,netFlux_Wm2)
copy_VarCoords(u_dTdx,u_dTdx_Wm2)
copy_VarCoords(v_dTdy,v_dTdy_Wm2)
copy_VarCoords(entrain,entrain_Wm2)

;====================================================================
; output the difference of some basic variables between two schemes
;====================================================================

T_grad_dx_diff		= T_grad_dx(1,:,:,:,:)-T_grad_dx(0,:,:,:,:)
copy_VarCoords_1(T_grad_dx(time|:,depth|:,lat|:,lon|:,case|:),T_grad_dx_diff)

T_grad_dy_diff		= T_grad_dy(1,:,:,:,:)-T_grad_dy(0,:,:,:,:)
copy_VarCoords_1(T_grad_dy(time|:,depth|:,lat|:,lon|:,case|:),T_grad_dy_diff)

Temp_mldmean_diff	= Temp_mldmean(1,:,:,:)-Temp_mldmean(0,:,:,:)
copy_VarCoords_1(Temp_mldmean(time|:,lat|:,lon|:,case|:),Temp_mldmean_diff)

u_mldmean_diff		= u_mldmean(1,:,:,:)-u_mldmean(0,:,:,:)
copy_VarCoords_1(u_mldmean(time|:,lat|:,lon|:,case|:),u_mldmean_diff)

v_mldmean_diff		= v_mldmean(1,:,:,:)-v_mldmean(0,:,:,:)
copy_VarCoords_1(v_mldmean(time|:,lat|:,lon|:,case|:),v_mldmean_diff)

w_mld_diff		= w_mld(1,:,:,:)-w_mld(0,:,:,:)
copy_VarCoords_1(w_mld(time|:,lat|:,lon|:,case|:),w_mld_diff)

deltaT_mld_diff		= deltaT_mld(1,:,:,:)-deltaT_mld(0,:,:,:)
copy_VarCoords_1(deltaT_mld(time|:,lat|:,lon|:,case|:),deltaT_mld_diff)

HMXL_diff		= HMXL(1,:,:,:)-HMXL(0,:,:,:)
copy_VarCoords_1(HMXL(time|:,lat|:,lon|:,case|:),HMXL_diff)

residual_Wm2_diff	= residual_Wm2(1,:,:,:)-residual_Wm2(0,:,:,:)
copy_VarCoords_1(residual_Wm2(time|:,lat|:,lon|:,case|:),residual_Wm2_diff)

sumall_Wm2_diff		= sumall_Wm2(1,:,:,:)-sumall_Wm2(0,:,:,:)
copy_VarCoords_1(sumall_Wm2(time|:,lat|:,lon|:,case|:),sumall_Wm2_diff)

netFlux_Wm2_diff	= netFlux_Wm2(1,:,:,:)-netFlux_Wm2(0,:,:,:)
copy_VarCoords_1(netFlux_Wm2(time|:,lat|:,lon|:,case|:),netFlux_Wm2_diff)

u_dTdx_Wm2_diff		= u_dTdx_Wm2(1,:,:,:)-u_dTdx_Wm2(0,:,:,:)
copy_VarCoords_1(u_dTdx_Wm2(time|:,lat|:,lon|:,case|:),u_dTdx_Wm2_diff)

v_dTdy_Wm2_diff		= v_dTdy_Wm2(1,:,:,:)-v_dTdy_Wm2(0,:,:,:)
copy_VarCoords_1(v_dTdy_Wm2(time|:,lat|:,lon|:,case|:),v_dTdy_Wm2_diff)

entrain_Wm2_diff	= entrain_Wm2(1,:,:,:)-entrain_Wm2(0,:,:,:)
copy_VarCoords_1(entrain_Wm2(time|:,lat|:,lon|:,case|:),entrain_Wm2_diff)


residual_diff		= residual(1,:,:,:)-residual(0,:,:,:)
copy_VarCoords_1(residual(time|:,lat|:,lon|:,case|:),residual_diff)

sumall_diff		= sumall(1,:,:,:)-sumall(0,:,:,:)
copy_VarCoords_1(sumall(time|:,lat|:,lon|:,case|:),sumall_diff)

netFlux_diff		= netFlux(1,:,:,:)-netFlux(0,:,:,:)
copy_VarCoords_1(netFlux(time|:,lat|:,lon|:,case|:),netFlux_diff)

u_dTdx_diff		= u_dTdx(1,:,:,:)-u_dTdx(0,:,:,:)
copy_VarCoords_1(u_dTdx(time|:,lat|:,lon|:,case|:),u_dTdx_diff)

v_dTdy_diff		= v_dTdy(1,:,:,:)-v_dTdy(0,:,:,:)
copy_VarCoords_1(v_dTdy(time|:,lat|:,lon|:,case|:),v_dTdy_diff)

entrain_diff		= entrain(1,:,:,:)-entrain(0,:,:,:)
copy_VarCoords_1(entrain(time|:,lat|:,lon|:,case|:),entrain_diff)

;========================================================
; write into one NC file to check them
;========================================================
system("/bin/rm -f ocn-heat-budget-evolution-terms"+"_"+season(is)+"-v1.nc")
ncdf = addfile("ocn-heat-budget-evolution-terms"+"_"+season(is)+"-v1.nc" ,"c")

ncdf->residual 		= residual
ncdf->sumall 		= sumall
ncdf->netFlux		= netFlux
ncdf->u_dTdx 		= u_dTdx
ncdf->v_dTdy 		= v_dTdy
ncdf->entrain		= entrain

ncdf->residual_diff 	= residual_diff
ncdf->sumall_diff 	= sumall_diff
ncdf->netFlux_diff	= netFlux_diff
ncdf->u_dTdx_diff	= u_dTdx_diff
ncdf->v_dTdy_diff	= v_dTdy_diff
ncdf->entrain_diff	= entrain_diff


;ncdf->residual_Wm2 	= residual_Wm2
;ncdf->sumall_Wm2 	= sumall_Wm2
;ncdf->netFlux_Wm2	= netFlux_Wm2
;ncdf->u_dTdx_Wm2 	= u_dTdx_Wm2
;ncdf->v_dTdy_Wm2 	= v_dTdy_Wm2
;ncdf->entrain_Wm2	= entrain_Wm2
;
;ncdf->residual_Wm2_diff = residual_Wm2_diff
;ncdf->sumall_Wm2_diff 	= sumall_Wm2_diff
;ncdf->netFlux_Wm2_diff	= netFlux_Wm2_diff
;ncdf->u_dTdx_Wm2_diff	= u_dTdx_Wm2_diff
;ncdf->v_dTdy_Wm2_diff	= v_dTdy_Wm2_diff
;ncdf->entrain_Wm2_diff	= entrain_Wm2_diff
;
;ncdf->finalvar_Wm2	= finalvar_Wm2
;ncdf->finalvardiff_Wm2	= finalvardiff_Wm2

system("/bin/rm -f ocn-heat-budget-evolution-not-terms"+"_"+season(is)+"-v1.nc")
ncdf2 = addfile("ocn-heat-budget-evolution-not-terms"+"_"+season(is)+"-v1.nc" ,"c")

;ncdf2->T_grad_dx		= T_grad_dx
;ncdf2->T_grad_dy		= T_grad_dy
ncdf2->Temp_mldmean		= Temp_mldmean
ncdf2->u_mldmean		= u_mldmean
ncdf2->v_mldmean		= v_mldmean
ncdf2->w_mld			= w_mld
ncdf2->deltaT_mld		= deltaT_mld
ncdf2->deltaT_mld_old		= deltaT_mld_old
ncdf2->HMXL			= HMXL

;ncdf->T_grad_dx_diff	= T_grad_dx_diff
;ncdf->T_grad_dy_diff	= T_grad_dy_diff
;ncdf->Temp_mldmean_diff	= Temp_mldmean_diff
;ncdf->u_mldmean_diff	= u_mldmean_diff
;ncdf->v_mldmean_diff	= v_mldmean_diff
;ncdf->w_mld_diff	= w_mld_diff
;ncdf->deltaT_mld_diff	= deltaT_mld_diff
;ncdf->HMXL_diff		= HMXL_diff
;
ncdf->index_mld		= index_mld

system("date")

end do; do is=
end do; do icompset=
end do; do icomp=
                               
end
