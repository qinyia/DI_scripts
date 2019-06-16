load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

load "./taylor_stats.ncl"

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ttest_all(xtmp[*][*]:numeric,ytmp[*][*]:numeric)
begin

xAve = dim_avg (xtmp)              ; calculate means at each grid point 
yAve = dim_avg (ytmp)
xVar = dim_variance (xtmp)         ; calculate variances
yVar = dim_variance (ytmp)

;print(xAve+"		"+yAve+"		"+xVar+"	"+yVar)

sigr = 0.05                        ; critical sig lvl for r

xN   = dimsizes(xtmp(0,:))
yN   = dimsizes(ytmp(0,:))

iflag= False                        ; population variance similar
prob = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False) 
prob@_FillValue		= -999.

printVarSummary(prob)

probf	= where(prob.lt.sigr,prob,prob@_FillValue)
print(probf)

return([/prob,probf/])

end


function regrid_diff_time(time[*]:numeric,\
                          var_def[*][*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		          var_new[*][*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

lat_spc			 	= fspan(-90,90,72)
lon_spc				= fspan(0,360,144)

;horizontal interpolation
var_interp 			= time_lat_lon_interpolate_2vars_specgrid (dimsizes(time), \
							 		   var_def,lat_def,lon_def, \
							 	 	   var_new,lat_new,lon_new, \
									   lat_spc,lon_spc)
var_interp_diff 		= var_interp(0,:,:,:)
var_interp_diff		 	= var_interp(1,:,:,:) - var_interp(0,:,:,:)

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin

case=(/"def","new"/)

;season=(/"ANN","DJF","MAM","JJA","SON"/)
season=(/"ANN"/)

var=(/"PRECT"/)

Dir				= "/home/lyl/WORK4/cesm1_2_1/archive/"

fdir_append 			= (/\
          		    "/atm/hist/", \
	  		    "/atm/hist/" /)

do ivar=0,dimsizes(var)-1


do icase=0,dimsizes(case)-1
print(case(icase))

if(case(icase).eq."def")then
	casename=(/"FC5_f19f19_MG15_cosp","BC5_f19g16_cosp"/)
	fdir=(/Dir,Dir/)
end if

if(case(icase).eq."new")then
	casename=(/"FC5_f19f19_MG15_mac3_cosp","BC5_f19g16_mac2_cosp"/)
	fdir=(/Dir,Dir/)
end if


;get_var_atm_TimeLevLatLon_addfiles(Dir:string,fdir_append:string,case:string,var:string,type:string)
; fils                    = systemfunc ("ls "+Dir+case+fdir_append+"*.nc") ; file paths

result_def			= get_var_atm_TimeLatLon_addfiles(fdir(0),fdir_append(0),casename(0),var(ivar),"model")
var_def				= result_def[0]
lat_def				= result_def[1]
lon_def				= result_def[2]
time_def			= result_def[3]

result_new			= get_var_atm_TimeLatLon_addfiles(fdir(1),fdir_append(1),casename(1),var(ivar),"model")
var_new				= result_new[0]
lat_new				= result_new[1]
lon_new				= result_new[2]

printVarSummary(var_def)
printVarSummary(var_new)

result 				= regrid_diff_time(time_def, \
			 			   var_def,lat_def,lon_def,\
			 			   var_new,lat_new,lon_new)

var_interp 			= result[0]
var_interp_diff 		= result[1]

dsizes				= dimsizes(var_interp)
dsizes_diff			= dimsizes(var_interp_diff)

delete([/var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/dsizes,dsizes_diff/])
;delete([/var_def_tmp,var_new_tmp/])

printVarSummary(var_interp) ; case, time, lat, lon

dsizes 				= dimsizes(var_interp)
ntime				= dsizes(1)

;-------------------------------------------
; calculate the PAI index
;-------------------------------------------

lat_here 			= var_interp&lat
rad 				= 4.0*atan(1.0)/180.
clat 				= cos(lat_here*rad)
clat!0 				= "lat"
clat&lat 			= lat_here

lonSa				= (/0,140,280,0/)
lonEa				= (/360,280,360,90/)

if(icase.eq.0)then
	var_interp_NH_all 		= new((/dimsizes(lonSa),4,ntime/),float) ; (lonSa, (amip_def,amip_new,cmip_def,cmip_new))
	var_interp_SH_all 		= new((/dimsizes(lonSa),4,ntime/),float)
	var_interp_tropical_all 	= new((/dimsizes(lonSa),4,ntime/),float)
	var_interp_PAI_all 		= new((/dimsizes(lonSa),4,ntime/),float)
end if

do ii=0,dimsizes(lonSa)-1
var_interp_NH 			= wgt_areaave_Wrap(var_interp(:,:,{0:20},{lonSa(ii):lonEa(ii)}),  clat({0:20}),1.0,1)
var_interp_SH 			= wgt_areaave_Wrap(var_interp(:,:,{-20:0},{lonSa(ii):lonEa(ii)}), clat({-20:0}),1.0,1)
var_interp_tropical 		= wgt_areaave_Wrap(var_interp(:,:,{-20:20},{lonSa(ii):lonEa(ii)}),clat({-20:20}),1.0,1)
var_interp_PAI 			= (var_interp_NH-var_interp_SH)/var_interp_tropical
;printVarSummary(var_interp_NH) ; (case,time)

;;; this varible is used to put all results from two cases.
var_interp_NH_all(ii,2*icase:2*icase+1,:) 		= var_interp_NH
var_interp_SH_all(ii,2*icase:2*icase+1,:) 		= var_interp_SH
var_interp_tropical_all(ii,2*icase:2*icase+1,:) 	= var_interp_tropical
var_interp_PAI_all(ii,2*icase:2*icase+1,:) 		= var_interp_PAI

delete(var_interp_NH)
delete(var_interp_SH)
delete(var_interp_tropical)
delete(var_interp_PAI)

end do ; do ii=
printVarSummary(var_interp_PAI_all)
end do; do icase=


xtmp			= var_interp_PAI_all(:,0,:);-conform_dims(dimsizes(var_interp_PAI_all(:,0,:)),dim_avg(var_interp_PAI_all(:,0,:)),0)
ytmp			= var_interp_PAI_all(:,2,:);-conform_dims(dimsizes(var_interp_PAI_all(:,2,:)),dim_avg(var_interp_PAI_all(:,2,:)),0)

;print(xtmp(0,:)+"	"+ytmp(0,:))
;exit

results			= ttest_all(xtmp,ytmp)
prob			= results[0]
probf			= results[1]
print(dim_avg(xtmp)+"	"+dim_avg(ytmp)+"		"+prob+"	"+probf)
exit


;-------------------------------------------
; calculate southern ITCZ (SI) index: 
; annual mean precipitation over (20S-0, 100-150W)
;-------------------------------------------
latS				= -20
latE				= 0
lonS				= 210
lonE				= 260
var_interp_SI 			= wgt_areaave(var_interp(:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1)

;;; this varible is used to put all results from two cases.
var_interp_SI_all 				= new((/4/),float)
var_interp_SI_all(2*icase:2*icase+1) 		= var_interp_SI
delete(var_interp_SI)

;-------------------------------------------
; calculate skill score defined by Taylor (2001)
; S = (1+R)^4/(4*(SDR+1/SDR)^2)
; R is a pattern correlation between the models and the observations,
; SDR is the ratio of spatial standard deviations of the models against that of observation.
; Thus, this score quantifies the similarity of the distribution and amplitude of a spatial pattern to the observations.
; 2018-05-26: add sectional skill score estimate
;-------------------------------------------
; R: pattern correlation
latSa				= (/-30,-30,-30,-30/)
latEa				= (/30,30,30,30/)
lonSa				= (/0,140,280,0/)
lonEa				= (/360,280,360,90/)

if(icase.eq.0)then
var_interp_Ra_all 		= new((/dimsizes(latSa),4/),float)
var_interp_ratio_all		= new((/dimsizes(latSa),4/),float)
end if

exit 
; I met problem here: there is no 144 samples for observation precipitation.
; so I cannot get the ratio of standard deviation and spatial correlation for two simulations.
; I am going to change to the monthly climatology.


do ii=0,dimsizes(latSa)-1

pr_amip				= var_interp(1,:,{latSa(ii):latEa(ii)},{lonSa(ii):lonEa(ii)})
pr_cmip				= var_interp(2,:,{latSa(ii):latEa(ii)},{lonSa(ii):lonEa(ii)})

; try to directly use the function from Taylor diagram: taylor_stats
;do it=0,ntime-1
tstVar				= pr_amip(it,:,:)
refVar				= pr_obs
stat_taylor			= taylor_stats(tstVar,refVar,clat({latSa(ii):latEa(ii)}),0)
Ra_amip				= stat_taylor(0) ;correlation
ratio_amip			= stat_taylor(1) ; ratio

tstVar				= pr_cmip
refVar				= pr_obs
stat_taylor			= taylor_stats(tstVar,refVar,clat({latSa(ii):latEa(ii)}),0)
Ra_cmip				= stat_taylor(0) ;correlation
ratio_cmip			= stat_taylor(1) ; ratio

;;; this varible is used to put all results from two cases together.
if(icase.eq.0)then
	var_interp_Ra_all(ii,3*icase:3*icase+2) 		= (/Ra_obs,Ra_amip,Ra_cmip/)
	var_interp_ratio_all(ii,3*icase:3*icase+2) 		= (/ratio_obs,ratio_amip,ratio_cmip/)
else
	var_interp_Ra_all(ii,3*icase:3*icase+1) 		= (/Ra_amip,Ra_cmip/)
	var_interp_ratio_all(ii,3*icase:3*icase+1) 		= (/ratio_amip,ratio_cmip/)
end if

delete([/pr_obs,pr_amip,pr_cmip/])
delete([/tstVar,refVar,stat_taylor,Ra_obs,ratio_obs,Ra_amip,ratio_amip,Ra_cmip,ratio_cmip/])
end do ; do ii=

; S = (1+R)^4/(4*(SDR+1/SDR)^2)
var_interp_S_all				= (1+var_interp_Ra_all)^4/(4.*(var_interp_ratio_all+1./var_interp_ratio_all)^2)

;end do; do icase=

printVarSummary(var_interp_S_all)

; -----------PAI index-----------------------------------------------
var_interp_PAI_all!0					= "region"
var_interp_PAI_all!1 					= "case"
var_interp_PAI_all&case 				= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)

; -----------SI index-----------------------------------------------
var_interp_SI_all!0 					= "case"
var_interp_SI_all&case 					= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)

; -----------skill score defined by Taylor (2001)-----------------------------------------------
var_interp_Ra_all!0				 	= "region"
var_interp_ratio_all!0				 	= "region"
var_interp_S_all!0				 	= "region"

var_interp_Ra_all!1 					= "case"
var_interp_Ra_all&case 					= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)
var_interp_ratio_all!1 					= "case"
var_interp_ratio_all&case 				= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)
var_interp_S_all!1 					= "case"
var_interp_S_all&case 					= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)

regions							= (/"global","Pacific(140E-80W)","South-America(80W-0)","Atlantic-Indian(0-90E)"/)

do ii=0,dimsizes(latSa)-1
print("--------------"+regions(ii)+"-----------------")
; ----------organize all indexes in one table------------------------------------------------
print("cases		PAI		SI		pat_cor		ratio		skill_score")
print(var_interp_Ra_all&case+"	"\
				+sprintf("%9.3f",var_interp_PAI_all(ii,:))+"	"\
				+sprintf("%9.3f",var_interp_SI_all)+"	"\
				+sprintf("%9.3f",var_interp_Ra_all(ii,:))+"	"\
				+sprintf("%9.3f",var_interp_ratio_all(ii,:))+"	"\
				+sprintf("%9.3f",var_interp_S_all(ii,:))+"	"\
     )

end do ; do ii=

end do ; do ivar=



                               
end
