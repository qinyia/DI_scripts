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
varTimeDepthLatLon	=(/"TEMP","UVEL","VVEL","WVEL","HOR_DIFF","VDC_T"/)

varTimeLatLon		=(/"HMXL","SHF","SHF_QSW"/)

compare_type 		= (/"cmip-obs"/)

case			= (/"def","new"/)

Output_TimeDepthLatLon		= True
Output_TimeLatLon		= True

factor			= 30.*24*3600.

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
depth!0			=	"depth"
depth&depth		=	depth

printVarSummary(HMXL)

; (vars,case,time,depth,lat,lon)
TEMP			= 	vars_TimeDepthLatLon(0,:,:,:,:,:)
UVEL			= 	vars_TimeDepthLatLon(1,:,:,:,:,:)
VVEL			= 	vars_TimeDepthLatLon(2,:,:,:,:,:)
WVEL			= 	vars_TimeDepthLatLon(3,:,:,:,:,:)
HOR_DIFF		= 	vars_TimeDepthLatLon(4,:,:,:,:,:) ; cm^2/s
HOR_DIFF		= 	HOR_DIFF*1e-4 ;m^2/s
VDC_T			= 	vars_TimeDepthLatLon(5,:,:,:,:,:) ; cm^2/s
VDC_T			=	VDC_T*1e-6	; m^2/s


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
deltaT_mld			= 	new((/ncase,ntime,nlat,nlon/),float)

; 2017-12-23
dTdz_mld			= 	new((/ncase,ntime,nlat,nlon/),float) ; K/m

VDC_T_mldmean			=	new((/ncase,ntime,nlat,nlon/),float) ; m^2/s
dTdz_mldmean			= 	new((/ncase,ntime,nlat,nlon/),float)
d2Tdz2_mldmean			= 	new((/ncase,ntime,nlat,nlon/),float)

; 2017-12-28
VDC_T_mld			=	new((/ncase,ntime,nlat,nlon/),float)
d2Tdz2_mld			=	new((/ncase,ntime,nlat,nlon/),float)

;========================================================
; calculate vertical advection
;========================================================
w_dTdz				=	new((/ncase,ntime,nlat,nlon/),float)

w_dTdz_lev	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)

; calculate vertical difference --> dTdz is generally negative because of the positive direction of depth is downward.
; so change its positive direction is upward to be consistent with other variables.
; WVEL is positive upward.
dTdz				= 	center_finite_diff_n(TEMP,depth,False,0,2)
dTdz				=	dTdz*(-1.)
copy_VarCoords(TEMP,dTdz)

; calculate vertical T Laplacian: d2Tdz2
d2Tdz2				=	center_finite_diff_n(dTdz,depth,False,0,2)
d2Tdz2				=	d2Tdz2*(-1.)
copy_VarCoords(TEMP,d2Tdz2)
printVarSummary(d2Tdz2)

; vertical advection
w_dTdz_lev			=	WVEL*dTdz
copy_VarCoords(WVEL,w_dTdz_lev)
printVarSummary(w_dTdz_lev)


;system("/bin/rm -f dTdz-v0.nc")
;ncdf 				=	addfile("dTdz-v0.nc" ,"c")
;
;ncdf->dTdz			= 	dTdz

;========================================================
; calculate horizontal advection: zonal and meridional && horizontal diffusion
;========================================================
; 2017-12-18 23:07:24 revised
; the calcualtion of u and v advection terms are based on each level and then get the vertical mean.
; so move the calculation of horizontal advection in front of the do loops.
u_dTdx				= 	new((/ncase,ntime,nlat,nlon/),float)
v_dTdy				= 	new((/ncase,ntime,nlat,nlon/),float)
;2017-12-23
hordiff				=	new((/ncase,ntime,nlat,nlon/),float)
vertdiff0			=	new((/ncase,ntime,nlat,nlon/),float)
vertdiff1			=	new((/ncase,ntime,nlat,nlon/),float) ; 2017-12-26
vertdiff2			=	new((/ncase,ntime,nlat,nlon/),float) ; 2017-12-28
vertdiff3			=	new((/ncase,ntime,nlat,nlon/),float) ; 2017-12-29
Qpenflux_mld			=	new((/ncase,ntime,nlat,nlon/),float) ; 2017-12-29

u_dTdx_lev	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
v_dTdy_lev	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
hordiff_lev			=	new((/ncase,ntime,ndepth,nlat,nlon/),float)
vertdiff1_lev			=	new((/ncase,ntime,ndepth,nlat,nlon/),float) ; 2017-12-26

T_grad_dx	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
T_grad_dy	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)

T_grad_dxdx	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
T_grad_dxdy	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
T_grad_dydx	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
T_grad_dydy	 		= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)


copy_VarCoords(TEMP,u_dTdx_lev)
copy_VarCoords(TEMP,v_dTdy_lev)
copy_VarCoords(TEMP,hordiff_lev)
copy_VarCoords(TEMP,vertdiff1_lev) ; 2017-12-26

copy_VarCoords(TEMP,T_grad_dx)
copy_VarCoords(TEMP,T_grad_dy)

printVarSummary(u_dTdx_lev)

;----------------------
; get latitudinal and longitudinal gradient
; grad_latlon_cfd only can deal with 2-4 dimensions variables. 
; so I have to do loop to calculate two cases' results.
;------------------------
; 2017-12-17: problem: T_grad_dx gets some strange values over north pole. needs to be solved.
; 2017-12-18: the problem comes from the regrided data after interpolation: fv1.9x2.5.
; I am still not sure why it has this problem. Now I change the grid into 1x1d before, and that problem disappears.

do icase=0,ncase-1

gradLatLon  			= 	grad_latlon_cfd (TEMP(icase,:,:,:,:), vars_TimeDepthLatLon&lat, vars_TimeDepthLatLon&lon, False, False)
T_grad_dy(icase,:,:,:,:)  	= 	gradLatLon[0] 
T_grad_dx(icase,:,:,:,:)  	= 	gradLatLon[1]
delete(gradLatLon)

;calculate Laplacian
gradLatLon  			= 	grad_latlon_cfd (T_grad_dx(icase,:,:,:,:), vars_TimeDepthLatLon&lat, vars_TimeDepthLatLon&lon, False, False)
T_grad_dxdy(icase,:,:,:,:)  	= 	gradLatLon[0] 
T_grad_dxdx(icase,:,:,:,:)  	= 	gradLatLon[1]
delete(gradLatLon)

gradLatLon  			= 	grad_latlon_cfd (T_grad_dy(icase,:,:,:,:), vars_TimeDepthLatLon&lat, vars_TimeDepthLatLon&lon, False, False)
T_grad_dydy(icase,:,:,:,:)  	= 	gradLatLon[0] 
T_grad_dydx(icase,:,:,:,:)  	= 	gradLatLon[1]
delete(gradLatLon)

printVarSummary(T_grad_dydy)


; zonal advection: u_dTdx
u_dTdx_lev(icase,:,:,:,:)	= 	T_grad_dx(icase,:,:,:,:)
u_dTdx_lev(icase,:,:,:,:)	= 	UVEL(icase,:,:,:,:)*T_grad_dx(icase,:,:,:,:)

; meridional advection: v_dTdy
v_dTdy_lev(icase,:,:,:,:)	= 	T_grad_dy(icase,:,:,:,:)
v_dTdy_lev(icase,:,:,:,:)	= 	VVEL(icase,:,:,:,:)*T_grad_dy(icase,:,:,:,:)

; horizontal diffusion
hordiff_lev(icase,:,:,:,:)	= 	T_grad_dxdx(icase,:,:,:,:)
hordiff_lev(icase,:,:,:,:)	= 	HOR_DIFF(icase,:,:,:,:)*(T_grad_dxdx(icase,:,:,:,:)+T_grad_dydy(icase,:,:,:,:))

end do ; do icase=

T_grad_dx@long_name 		= 	"zonal gradient"
T_grad_dy@long_name 		= 	"meridional gradient"
T_grad_dx@units     		= 	"degC/m"
T_grad_dy@units     		= 	"degC/m"

printVarSummary(u_dTdx_lev)
printVarSummary(T_grad_dx)
printVarSummary(hordiff_lev)
;========================================================
; vertical diffusion_v1, by d2Tdz2*VDC_T
;========================================================
vertdiff1_lev			=	d2Tdz2*VDC_T
vertdiff1_lev			= 	vertdiff1_lev*factor
copy_VarCoords(VDC_T,vertdiff1_lev)
printVarSummary(vertdiff1_lev)
;========================================================
; calculate penetrated SW flux
; estimated by using the empirical formula of Paulson and Simpson (1977)
; referred from Liu et al., 2012
; parameters are following the Rosati and Miyakoda (1988)
;========================================================
R					=	0.58
z1					=	0.35 ; m
z2					=	23 ; m
SHF_QSW					=	vars_TimeLatLon(2,:,:,:,:)	
Qpenflux_lev				= 	new((/ncase,ntime,ndepth,nlat,nlon/),float)
copy_VarCoords(TEMP,Qpenflux_lev)
do ilev=0,ndepth-1
Qpenflux_lev(:,:,ilev,:,:)		= 	SHF_QSW*(R*exp(-1.*depth(ilev)/z1)+(1-R)*exp(-1.*depth(ilev)/z2))
end do

system("/bin/rm -f d2Tdz2-v0.nc")
ncdf 				=	addfile("d2Tdz2-v0.nc" ,"c")

ncdf->u_dTdx_lev		=	u_dTdx_lev
ncdf->v_dTdy_lev		=	v_dTdy_lev

ncdf->d2Tdz2			=	d2Tdz2
ncdf->VDC_T			=	VDC_T
ncdf->dTdz			=	dTdz

ncdf->Qpenflux_lev		=	Qpenflux_lev

;========================================================
; Get Mixed-layer mean of all necessary variables
;========================================================

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
;do itime = 0,0
do ilat = 0,nlat-1
do ilon = 0,nlon-1


index_mld_tmp(icase,itime,ilat,ilon) 			= 	closest_val(HMXL(icase,itime,ilat,ilon),depth)	

index_mld(icase,itime,ilat,ilon)			=	index_mld_tmp(icase,itime,ilat,ilon)

; 2017-12-25: be careful about this choice: should be "ge" not "gt". "gt" will lose the situation with index=1
;if(.not.ismissing(index_mld_tmp(icase,itime,ilat,ilon)).and.index_mld_tmp(icase,itime,ilat,ilon).ge.1)then
;	if(depth(index_mld_tmp(icase,itime,ilat,ilon)).gt.hmxl(icase,itime,ilat,ilon))then
;		index_mld(icase,itime,ilat,ilon)	=	index_mld_tmp(icase,itime,ilat,ilon)-1
;	else
;		index_mld(icase,itime,ilat,ilon)	=	index_mld_tmp(icase,itime,ilat,ilon)
;	end if
;end if

print(icase+"  "+itime+"  "+ilat+"  "+ilon+" index="+index_mld(icase,itime,ilat,ilon))

if(.not.ismissing(index_mld(icase,itime,ilat,ilon)))then
	; vertical mean of all inputted variables
	vars_TimeDepthLatLon_mldmean(:,icase,itime,ilat,ilon) 	= dim_avg_wgt_n_Wrap(vars_TimeDepthLatLon(:,icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,1)

	; vertical mean of dTdz
	dTdz_mldmean(icase,itime,ilat,ilon)			= dim_avg_wgt_n_Wrap(dTdz(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	; vertical mean of d2Tdz2
	d2Tdz2_mldmean(icase,itime,ilat,ilon)			= dim_avg_wgt_n_Wrap(d2Tdz2(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	; vertical velocity at mixed-layer depth
	w_mld(icase,itime,ilat,ilon)				= (/vars_TimeDepthLatLon(3,icase,itime,index_mld(icase,itime,ilat,ilon),ilat,ilon)/)
	; vertical temperature gradient at mixed-layer depth
	deltaT_mld(icase,itime,ilat,ilon)			= (/vars_TimeDepthLatLon_mldmean(0,icase,itime,ilat,ilon)-vars_TimeDepthLatLon(0,icase,itime,index_mld(icase,itime,ilat,ilon)+1,ilat,ilon)/)

	; vertical mean of horizontal advection
	; 2017-12-19 revised
	u_dTdx(icase,itime,ilat,ilon)				= dim_avg_wgt_n_Wrap(u_dTdx_lev(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	v_dTdy(icase,itime,ilat,ilon)				= dim_avg_wgt_n_Wrap(v_dTdy_lev(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	; vertical mean of vertical advection
	w_dTdz(icase,itime,ilat,ilon)				= dim_avg_wgt_n_Wrap(w_dTdz_lev(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	; vertical mean of horizontal diffusion
	hordiff(icase,itime,ilat,ilon)				= dim_avg_wgt_n_Wrap(hordiff_lev(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	; vertical mean of vertical diffusion_v1
	vertdiff1(icase,itime,ilat,ilon)			= dim_avg_wgt_n_Wrap(vertdiff1_lev(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	; dTdz at mixed layer base
;	dTdz_mld(icase,itime,ilat,ilon)			 	= (TEMP(icase,itime,index_mld(icase,itime,ilat,ilon),ilat,ilon)-TEMP(icase,itime,index_mld(icase,itime,ilat,ilon)+1,ilat,ilon))/(depth(index_mld(icase,itime,ilat,ilon))-depth(index_mld(icase,itime,ilat,ilon)+1))
	dTdz_mld(icase,itime,ilat,ilon)				= dTdz(icase,itime,index_mld(icase,itime,ilat,ilon),ilat,ilon)
	print("dTdz_mld="+dTdz_mld(icase,itime,ilat,ilon)+"  dTdz_mldmean="+dTdz_mldmean(icase,itime,ilat,ilon))
	; vertical mean of vertical diffusivity coefficient
	VDC_T_mldmean(icase,itime,ilat,ilon)			= dim_avg_wgt_n_Wrap(VDC_T(icase,itime,:index_mld(icase,itime,ilat,ilon),ilat,ilon),wgt(:index_mld(icase,itime,ilat,ilon)),1,0)
	; d2Tdz2 at mixed layer base
	d2Tdz2_mld(icase,itime,ilat,ilon)			= d2Tdz2(icase,itime,index_mld(icase,itime,ilat,ilon),ilat,ilon)
	; VDC_T at mixed layer base
	VDC_T_mld(icase,itime,ilat,ilon)			= VDC_T(icase,itime,index_mld(icase,itime,ilat,ilon),ilat,ilon)
	; Qpenflux_lev at mixed layer base
	Qpenflux_mld(icase,itime,ilat,ilon)				= Qpenflux_lev(icase,itime,index_mld(icase,itime,ilat,ilon),ilat,ilon)
else
	vars_TimeDepthLatLon_mldmean(:,icase,itime,ilat,ilon) 	=	vars_TimeDepthLatLon@_FillValue
	w_mld(icase,itime,ilat,ilon)				= 	vars_TimeDepthLatLon@_FillValue
	deltaT_mld(icase,itime,ilat,ilon)			=	vars_TimeDepthLatLon@_FillValue
	u_dTdx(icase,itime,ilat,ilon)				=	vars_TimeDepthLatLon@_FillValue
	v_dTdy(icase,itime,ilat,ilon)				=	vars_TimeDepthLatLon@_FillValue
	hordiff(icase,itime,ilat,ilon)				= 	vars_TimeDepthLatLon@_FillValue
	vertdiff1(icase,itime,ilat,ilon)			= 	vars_TimeDepthLatLon@_FillValue
	dTdz_mld(icase,itime,ilat,ilon)				= 	vars_TimeDepthLatLon@_FillValue
	VDC_T_mldmean(icase,itime,ilat,ilon)			=	vars_TimeDepthLatLon@_FillValue
	dTdz_mldmean(icase,itime,ilat,ilon)			=	vars_TimeDepthLatLon@_FillValue
	w_dTdz(icase,itime,ilat,ilon)				=	vars_TimeDepthLatLon@_FillValue
	d2Tdz2_mld(icase,itime,ilat,ilon)			=	vars_TimeDepthLatLon@_FillValue
	VDC_T_mld(icase,itime,ilat,ilon)			=	vars_TimeDepthLatLon@_FillValue
	Qpenflux_mld(icase,itime,ilat,ilon)			=	vars_TimeDepthLatLon@_FillValue
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
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),deltaT_mld)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),index_mld)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),u_dTdx)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),v_dTdy)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),hordiff)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),dTdz_mld)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),VDC_T_mldmean)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),dTdz_mldmean)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),w_dTdz)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),vertdiff1)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),d2Tdz2_mld)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),VDC_T_mld)
copy_VarCoords_1(vars_TimeDepthLatLon_mldmean(case|:,time|:,lat|:,lon|:,vars|:),Qpenflux_mld)

printVarSummary(vars_TimeDepthLatLon_mldmean)

Temp_mldmean 		= 	vars_TimeDepthLatLon_mldmean(0,:,:,:,:)
u_mldmean		= 	vars_TimeDepthLatLon_mldmean(1,:,:,:,:)
v_mldmean		= 	vars_TimeDepthLatLon_mldmean(2,:,:,:,:)

;========================================================
; calculate entrainment term
;========================================================

entrain			= HMXL
entrain 		= w_mld*deltaT_mld/HMXL
entrain			= entrain*(-1.0)

entrain 		= entrain*factor	;change from K/s to K/month

entrain@long_name	= "entrainment term"	
entrain@units		= "degC/month"

printVarSummary(entrain)

;========================================================
; calculate horizontal diffusion
;========================================================
hordiff			=	hordiff*factor

hordiff@long_name	=	"horizontal diffusion term: hordiff"
hordiff@units		=	"degC/month"

;========================================================
; calculate vertical diffusion
;========================================================
vertdiff0		= 	HMXL
vertdiff0		=	0.
vertdiff0		= 	VDC_T_mldmean*dTdz_mld/HMXL
vertdiff0		=	vertdiff0*factor*(-1.)

vertdiff0@long_name	=	"vertical diffusion term: vertdiff_v0 by divided mixed layer depth"
vertdiff0@units		=	"degC/month"
; pay attention: it seems the vertdiff is larger than other terms. need to double check.
printVarSummary(vertdiff0)
;----------------------------------------------------------------
;vertdiff1		=	VDC_T_mldmean*d2Tdz2_mldmean
vertdiff1		=	vertdiff1*factor*(-1.)

vertdiff1@long_name	=	"vertical diffusion term: vertically integrated VDC_T*d2Tdz2"
vertdiff1@units		=	"degC/month"
;----------------------------------------------------------------
vertdiff2		=	VDC_T_mld*d2Tdz2_mld
vertdiff2		=	vertdiff2*factor*(-1.)

vertdiff2@long_name	=	"vertical diffusion term:VDC_T_mld*d2Tdz2_mld" 
vertdiff2@units		=	"degC/month"
;----------------------------------------------------------------
vertdiff3		=	HMXL
vertdiff3		=	VDC_T_mld*dTdz_mld/HMXL
vertdiff3		=	vertdiff3*factor*(-1.)

vertdiff3@long_name	=	"vertical diffusion term: VDC_T_mld*dTdz_mld/HMXL"
vertdiff3@units		=	"degC/month"
;========================================================
; calculate horizontal advection: zonal and meridional
;========================================================
u_dTdx			=	u_dTdx*(-1.)
;change unit from K/s to K/month
u_dTdx 			= 	u_dTdx*factor

u_dTdx@long_name	= 	"advection term: u_dTdx"
u_dTdx@units 		= 	"degC/month"


; meridional advection: v_dTdy
v_dTdy			= 	v_dTdy*(-1.)
;change units
v_dTdy 			= 	v_dTdy*factor

v_dTdy@longname 	= 	"advection term: v_dTdy"
v_dTdy@units 		= 	"degC/month"

printVarSummary(hordiff)
;========================================================
; calculate vertical advection
;========================================================
w_dTdz			=	w_dTdz*factor*(-1.)

w_dTdz@longname		=	"vertical advection:w_dTdz"
w_dTdz@units		=	"degC/month"

;========================================================
; calculate the net energy flux term
;========================================================

rho_sw 			= 1030.		; kg/m3
cp_sw			= 4.0e3		; J/kg/K

SHF			= vars_TimeLatLon(1,:,:,:,:) ; total surface heat flux: (case,lat,lon)
netFlux 		= SHF/rho_sw/cp_sw/HMXL
copy_VarCoords(SHF,netFlux)

netFlux			= netFlux*factor
netFlux@longname 	= "net energy flux term"
netFlux@units		= "degC/month"

printVarSummary(netFlux)

;========================================================
; calculate the penetrative SW flux
;========================================================
Qpenflux		= Qpenflux_mld
Qpenflux		= Qpenflux_mld/rho_sw/cp_sw/HMXL
Qpenflux		= Qpenflux*factor*(-1.)
Qpenflux@longname	= "penetrative SW flux"
Qpenflux@units	= "degC/month"

printVarSummary(Qpenflux)

;========================================================
; put all terms into one array to facilitate the plots
;========================================================
Terms		= (/"netFlux","u_dTdx","v_dTdy","w_dTdz","entrain","hordiff","vertdiff0","vertdiff1","vertdiff2","vertdiff3","Qpenflux"/)

TermsDiff	= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","w_dTdz_diff","entrain_diff","hordiff_diff","vertdiff0_diff","vertdiff1_diff","vertdiff2_diff","vertdiff3_diff","Qpenflux_diff"/)

terms			= ispan(0,dimsizes(Terms)-1,1)
finalvar		= new((/dimsizes(Terms),ncase,ntime,nlat,nlon/),float)
finalvar!0		= "terms"
finalvar&terms		= terms

finalvar(0,:,:,:,:)	= netFlux
finalvar(1:,:,:,:,:)	= (/u_dTdx,v_dTdy,w_dTdz,entrain,hordiff,vertdiff0,vertdiff1,vertdiff2,vertdiff3,Qpenflux/)


finalvardiff 		= new((/dimsizes(Terms),ntime,nlat,nlon/),float)
finalvardiff		= finalvar(:,1,:,:,:)-finalvar(:,0,:,:,:)
copy_VarCoords_1(finalvar(terms|:,time|:,lat|:,lon|:,case|:),finalvardiff)

printVarSummary(finalvar)
printVarSummary(finalvardiff)
printMinMax(finalvar,True)
printMinMax(finalvardiff,True)

;========================================================
; write into one NC file to check them
;========================================================
system("/bin/rm -f test-cal-terms-vertdiff3-Qpenflux.nc")
ncdf = addfile("test-cal-terms-vertdiff3-Qpenflux.nc" ,"c")

nterms			= dimsizes(Terms)

do iterm=0,nterms-1
	ncdf->$Terms(iterm)$		=	finalvar(iterm,:,:,:,:)
end do

do iterm=0,nterms-1
	ncdf->$TermsDiff(iterm)$	=	finalvardiff(iterm,:,:,:)
end do


system("/bin/rm -f test-cal-not-terms-vertdiff3-Qpenflux.nc")
ncdf2 = addfile("test-cal-not-terms-vertdiff3-Qpenflux.nc" ,"c")

ncdf2->Temp_mldmean		= Temp_mldmean
ncdf2->u_mldmean		= u_mldmean
ncdf2->v_mldmean		= v_mldmean
ncdf2->w_mld			= w_mld
ncdf2->deltaT_mld		= deltaT_mld
ncdf2->HMXL			= HMXL
ncdf2->VDC_T_mldmean		= VDC_T_mldmean
ncdf2->dTdz_mld			= dTdz_mld

ncdf->index_mld			= index_mld

system("date")

end do; do is=
end do; do icompset=
end do; do icomp=
                               
end
