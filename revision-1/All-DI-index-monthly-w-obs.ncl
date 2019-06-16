load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

load "./taylor_stats.ncl"

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ttest_2D(xtmp[*][*]:numeric,ytmp[*][*]:numeric)
begin

xAve = dim_avg (xtmp)              ; calculate means at each grid point 
yAve = dim_avg (ytmp)
xVar = dim_variance (xtmp)         ; calculate variances
yVar = dim_variance (ytmp)

sigr = 0.05                        ; critical sig lvl for r

xN   = dimsizes(xtmp(0,:))
yN   = dimsizes(ytmp(0,:))

iflag= False                        ; population variance similar
prob = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False) 
prob@_FillValue		= -999.

printVarSummary(prob)

probf	= where(prob.lt.sigr,prob,prob@_FillValue)
print(probf)

delete([/xAve,yAve,xVar,yVar,sigr,xN,yN/])
return([/prob,probf/])

end

function ttest_1D(xtmp[*]:numeric,ytmp[*]:numeric)
begin

xAve = dim_avg (xtmp)              ; calculate means at each grid point 
yAve = dim_avg (ytmp)
xVar = dim_variance (xtmp)         ; calculate variances
yVar = dim_variance (ytmp)

sigr = 0.05                        ; critical sig lvl for r

xN   = dimsizes(xtmp)
yN   = dimsizes(ytmp)

;xN   	= 144
;yN	= 144

iflag= False                        ; population variance similar
prob = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False) 
prob@_FillValue		= -999.

printVarSummary(prob)

probf	= where(prob.lt.sigr,prob,prob@_FillValue)
print(probf)

delete([/xAve,yAve,xVar,yVar,sigr,xN,yN/])
return([/prob,probf/])

end



;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff_time(time[*]:numeric,\
		     var_obs[*][*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                     var_amip[*][*][*]:numeric,lat_amip[*]:numeric,lon_amip[*]:numeric,\
		     var_cmip[*][*][*]:numeric,lat_cmip[*]:numeric,lon_cmip[*]:numeric)
begin

lat_spc			 	= fspan(-90,90,72)
lon_spc				= fspan(0,360,144)


;horizontal interpolation
var_interp 			= time_lat_lon_interpolate_3vars_specgrid (dimsizes(time), \
							 		   var_obs,lat_obs,lon_obs, \
							 		   var_amip,lat_amip,lon_amip, \
							 	 	   var_cmip,lat_cmip,lon_cmip, \
									   lat_spc,lon_spc)
var_interp_diff = var_interp
var_interp_diff(0,:,:,:) = var_interp(1,:,:,:) - var_interp(0,:,:,:) ; def-obs
var_interp_diff(1,:,:,:) = var_interp(2,:,:,:) - var_interp(0,:,:,:) ; new-obs
var_interp_diff(2,:,:,:) = var_interp(2,:,:,:) - var_interp(1,:,:,:) ; new-obs

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



do ivar=0,dimsizes(var)-1

do icase=0,dimsizes(case)-1
print(case(icase))

if(case(icase).eq."def")then
	casename=(/"GPCP","FC5_f19f19_MG15_cosp","BC5_f19g16_cosp"/)
	fdir=(/"/home/share3/lyl/work3/qinyi/double-ITCZ/data/obs/GPCP_extend_12yr/",\
	       "/home/lyl/WORK4/cesm1_2_1/archive/",\
	       "/home/lyl/WORK4/cesm1_2_1/archive/"/)
end if

if(case(icase).eq."new")then
	casename=(/"GPCP","FC5_f19f19_MG15_mac3_cosp","BC5_f19g16_mac2_cosp"/)
	fdir=(/"/home/share3/lyl/work3/qinyi/double-ITCZ/data/obs/GPCP_extend_12yr/",\
	       "/home/lyl/WORK4/cesm1_2_1/archive/",\
	       "/home/lyl/WORK4/cesm1_2_1/archive/"/)
end if

fdir_append 			= (/"",\
          		    "/atm/hist/", \
	  		    "/atm/hist/" /)

result_obs 			= get_var_atm_TimeLatLon_addfiles(fdir(0),fdir_append(0),casename(0),var(ivar),"obs")
var_obs 			= result_obs[0]
lat_obs   			= result_obs[1]
lon_obs   			= result_obs[2]

result_amip 			= get_var_atm_TimeLatLon_addfiles(fdir(1),fdir_append(1),casename(1),var(ivar),"model")
var_amip 			= result_amip[0]
lat_amip   			= result_amip[1]
lon_amip   			= result_amip[2]
time				= result_amip[3]

result_cmip 			= get_var_atm_TimeLatLon_addfiles(fdir(2),fdir_append(2),casename(2),var(ivar),"model")
var_cmip 			= result_cmip[0]
lat_cmip   			= result_cmip[1]
lon_cmip   			= result_cmip[2]


result 				= regrid_diff_time(time,\
					      var_obs,lat_obs,lon_obs,\
	             			      var_amip,lat_amip,lon_amip,\
		     			      var_cmip,lat_cmip,lon_cmip)
var_interp 			= result[0] ; obs,amip,cmip
var_interp_diff 		= result[1] ; obs,amip,cmip

printVarSummary(var_interp)			; (case,time,lat,lon) = (3, 12, 72, 144)

dsizes				= dimsizes(var_interp)
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
	var_interp_NH_all 				= new((/dimsizes(lonSa),5,ntime/),float)
	var_interp_SH_all 				= new((/dimsizes(lonSa),5,ntime/),float)
	var_interp_tropical_all 			= new((/dimsizes(lonSa),5,ntime/),float)
	var_interp_PAI_all 				= new((/dimsizes(lonSa),5,ntime/),float)
end if

do ii=0,dimsizes(lonSa)-1
var_interp_NH 			= wgt_areaave_Wrap(var_interp(:,:,{0:20},{lonSa(ii):lonEa(ii)}),  clat({0:20}),1.0,1)
var_interp_SH 			= wgt_areaave_Wrap(var_interp(:,:,{-20:0},{lonSa(ii):lonEa(ii)}), clat({-20:0}),1.0,1)
var_interp_tropical 		= wgt_areaave_Wrap(var_interp(:,:,{-20:20},{lonSa(ii):lonEa(ii)}),clat({-20:20}),1.0,1)
var_interp_PAI 			= (var_interp_NH-var_interp_SH)/var_interp_tropical

;;; this varible is used to put all results from two cases.
if(icase.eq.0)then
	var_interp_NH_all(ii,3*icase:3*icase+2,:) 		= var_interp_NH
	var_interp_SH_all(ii,3*icase:3*icase+2,:) 		= var_interp_SH
	var_interp_tropical_all(ii,3*icase:3*icase+2,:) 	= var_interp_tropical
	var_interp_PAI_all(ii,3*icase:3*icase+2,:) 		= var_interp_PAI
else
	var_interp_NH_all(ii,3*icase:3*icase+1,:) 		= var_interp_NH(1:2,:)
	var_interp_SH_all(ii,3*icase:3*icase+1,:) 		= var_interp_SH(1:2,:)
	var_interp_tropical_all(ii,3*icase:3*icase+1,:) 	= var_interp_tropical(1:2,:)
	var_interp_PAI_all(ii,3*icase:3*icase+1,:) 		= var_interp_PAI(1:2,:)
end if
delete(var_interp_NH)
delete(var_interp_SH)
delete(var_interp_tropical)
delete(var_interp_PAI)

end do ; do ii=

;var_interp_PAI_all!0		 	= "region"
;var_interp_PAI_all&region		= fspan(0,dimsizes(lonSa)-1,dimsizes(lonSa))
;var_interp_PAI_all!1			= "case"
;var_interp_PAI_all&case			= fspan(0,4,5)
;var_interp_PAI_all!2			= "time"
;var_interp_PAI_all&time			= fspan(0,ntime-1,ntime)


printVarSummary(var_interp_PAI_all)
print(var_interp_PAI_all(0,0,:))

wks			= gsn_open_wks("x11","test")

res			= True
res@xyLineColors	= (/"black","red","blue"/)
res@xyLineThicknessF	= 3.0

plot			= gsn_y(wks,var_interp_PAI_all(0,:,:),res)
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
if(icase.eq.0)then
	var_interp_SI_all 				= new((/5,ntime/),float)
	var_interp_SI_all(3*icase:3*icase+2,:) 		= var_interp_SI
else
	var_interp_SI_all(3*icase:3*icase+1,:) 		= var_interp_SI(1:2,:)
end if
delete(var_interp_SI)

var_interp_SI_all!0			= "case"
var_interp_SI_all&case			= fspan(0,4,5)
var_interp_SI_all!1			= "time"
var_interp_SI_all&time			= fspan(0,ntime-1,ntime)

printVarSummary(var_interp_SI_all)
print(var_interp_SI_all(0,:))

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
var_interp_Ra_all 		= new((/dimsizes(latSa),5,ntime/),float)
var_interp_ratio_all		= new((/dimsizes(latSa),5,ntime/),float)
end if


do ii=0,dimsizes(latSa)-1

pr_obs				= var_interp(0,:,{latSa(ii):latEa(ii)},{lonSa(ii):lonEa(ii)})
pr_amip				= var_interp(1,:,{latSa(ii):latEa(ii)},{lonSa(ii):lonEa(ii)})
pr_cmip				= var_interp(2,:,{latSa(ii):latEa(ii)},{lonSa(ii):lonEa(ii)})

; try to directly use the function from Taylor diagram: taylor_stats
do it = 0,ntime-1
tstVar				= pr_obs(it,:,:)
refVar				= pr_obs(it,:,:)
stat_taylor			= taylor_stats(tstVar,refVar,clat({latSa(ii):latEa(ii)}),0)
Ra_obs				= stat_taylor(0) ;correlation
ratio_obs			= stat_taylor(1) ; ratio

tstVar				= pr_amip(it,:,:)
refVar				= pr_obs(it,:,:)
stat_taylor			= taylor_stats(tstVar,refVar,clat({latSa(ii):latEa(ii)}),0)
Ra_amip				= stat_taylor(0) ;correlation
ratio_amip			= stat_taylor(1) ; ratio

tstVar				= pr_cmip(it,:,:)
refVar				= pr_obs(it,:,:)
stat_taylor			= taylor_stats(tstVar,refVar,clat({latSa(ii):latEa(ii)}),0)
Ra_cmip				= stat_taylor(0) ;correlation
ratio_cmip			= stat_taylor(1) ; ratio

;;; this varible is used to put all results from two cases together.
if(icase.eq.0)then
	var_interp_Ra_all(ii,3*icase:3*icase+2,it) 		= (/Ra_obs,Ra_amip,Ra_cmip/)
	var_interp_ratio_all(ii,3*icase:3*icase+2,it) 		= (/ratio_obs,ratio_amip,ratio_cmip/)
else
	var_interp_Ra_all(ii,3*icase:3*icase+1,it) 		= (/Ra_amip,Ra_cmip/)
	var_interp_ratio_all(ii,3*icase:3*icase+1,it) 		= (/ratio_amip,ratio_cmip/)
end if

end do ; do it=

delete([/pr_obs,pr_amip,pr_cmip/])
delete([/tstVar,refVar,stat_taylor,Ra_obs,ratio_obs,Ra_amip,ratio_amip,Ra_cmip,ratio_cmip/])
end do ; do ii=

; S = (1+R)^4/(4*(SDR+1/SDR)^2)
var_interp_SS_all				= (1+var_interp_Ra_all)^4/(4.*(var_interp_ratio_all+1./var_interp_ratio_all)^2)

var_interp_SS_all!0		 	= "region"
var_interp_SS_all&region			= fspan(0,dimsizes(lonSa)-1,dimsizes(lonSa))
var_interp_SS_all!1			= "case"
var_interp_SS_all&case			= fspan(0,4,5)
var_interp_SS_all!2			= "time"
var_interp_SS_all&time			= fspan(0,ntime-1,ntime)

end do; do icase=
printVarSummary(var_interp_SS_all)
print(var_interp_SS_all(0,1,:))

;+++++++++++++++++++++++++++++++++++++++++ T test ++++++++++++++++++++++++++++++
; ttest to identify whether two simulations show very significant difference.

; ======================== PAI index ===================
xydiff_PAI			= new((/dimsizes(lonSa),6/),float)
xydiff_PAI!0			= "region"
xydiff_PAI&region		= (/"GLB","PA","SA","AI"/)
xydiff_PAI!1			= "case"
xydiff_PAI&case			= (/"amip-def","amip-new","sig-amip","cmip-def","cmip-new","sig-cmip"/)

do icomp=0,1 ; =0: amip; =1: cmip
xdiff_PAI			= var_interp_PAI_all(:,icomp+1,:)-var_interp_PAI_all(:,0,:)
ydiff_PAI			= var_interp_PAI_all(:,icomp+3,:)-var_interp_PAI_all(:,0,:)
results				= ttest_2D(xdiff_PAI,ydiff_PAI)
prob_PAI			= results[0]
probf_PAI			= results[1]
print(dim_avg(xdiff_PAI)+"	"+dim_avg(ydiff_PAI)+"		"+prob_PAI+"	"+probf_PAI)

xydiff_PAI(:,3*icomp)		= dim_avg(xdiff_PAI)
xydiff_PAI(:,3*icomp+1)		= dim_avg(ydiff_PAI)
xydiff_PAI(:,3*icomp+2)		= prob_PAI

end do ; do icomp=
print(xydiff_PAI(0,:))


; ======================== SI index ===================

xydiff_SI			= new((/6/),float)
xydiff_SI!0			= "case"
xydiff_SI&case			= (/"amip-def","amip-new","sig-amip","cmip-def","cmip-new","sig-cmip"/)

do icomp=0,1 ; =0: amip; =1: cmip
xdiff_SI			= var_interp_SI_all(icomp+1,:)-var_interp_SI_all(0,:)
ydiff_SI			= var_interp_SI_all(icomp+3,:)-var_interp_SI_all(0,:)
results				= ttest_1D(xdiff_SI,ydiff_SI)
prob_SI				= results[0]
probf_SI			= results[1]
print(dim_avg(xdiff_SI)+"	"+dim_avg(ydiff_SI)+"		"+prob_SI+"	"+probf_SI)

xydiff_SI(3*icomp)		= dim_avg(xdiff_SI)
xydiff_SI(3*icomp+1)		= dim_avg(ydiff_SI)
xydiff_SI(3*icomp+2)		= prob_SI

end do ; do icomp=
print(xydiff_SI)

; ======================== Skill score index ===================

xydiff_SS			= new((/dimsizes(lonSa),6/),float)
xydiff_SS!0			= "region"
xydiff_SS&region		= (/"GLB","PA","SA","AI"/)
xydiff_SS!1			= "case"
xydiff_SS&case			= (/"amip-def","amip-new","sig-amip","cmip-def","cmip-new","sig-cmip"/)

do icomp=0,1 ; =0: amip; =1: cmip
xdiff_SS			= var_interp_SS_all(:,icomp+1,:)-var_interp_SS_all(:,0,:)
ydiff_SS			= var_interp_SS_all(:,icomp+3,:)-var_interp_SS_all(:,0,:)
results				= ttest_2D(xdiff_SS,ydiff_SS)
prob_SS				= results[0]
probf_SS			= results[1]
print(dim_avg(xdiff_SS)+"	"+dim_avg(ydiff_SS)+"		"+prob_SS+"	"+probf_SS)

xydiff_SS(:,3*icomp)		= dim_avg(xdiff_SS)
xydiff_SS(:,3*icomp+1)		= dim_avg(ydiff_SS)
xydiff_SS(:,3*icomp+2)		= prob_SS

end do ; do icomp=
print(xydiff_SS(0,:))


regions							= (/"global","Pacific(140E-80W)","South-America(80W-0)","Atlantic-Indian(0-90E)"/)
;

do ii=0,dimsizes(latSa)-1

print("--------------"+regions(ii)+"-----------------")
; ----------organize all indexes in one table------------------------------------------------
print("cases,		PAI,		SI,		skill_score,")
print(xydiff_PAI&case+",	"\
				+sprintf("%9.5f",xydiff_PAI(ii,:))+",	"\
				+sprintf("%9.5f",xydiff_SI)+",	"\
				+sprintf("%9.5f",xydiff_SS(ii,:))+",	"\
     )

var_interp_PAI_all!1					= "case"
var_interp_PAI_all&case 				= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)

print("--------------"+regions(ii)+"-----------------")
; ----------organize all indexes in one table------------------------------------------------
print("cases,		PAI,		SI,		pat_cor,		ratio,		skill_score,")
print(var_interp_PAI_all&case+",	"\
				+sprintf("%9.3f",dim_avg(var_interp_PAI_all(ii,:,:)))+",	"\
				+sprintf("%9.3f",dim_avg(var_interp_SI_all(:,:)))+",	"\
				+sprintf("%9.3f",dim_avg(var_interp_Ra_all(ii,:,:)))+",	"\
				+sprintf("%9.3f",dim_avg(var_interp_ratio_all(ii,:,:)))+",	"\
				+sprintf("%9.3f",dim_avg(var_interp_SS_all(ii,:,:)))+",	"\
     )


end do ; do ii=



;; -----------PAI index-----------------------------------------------
;var_interp_PAI_all!0					= "region"
;var_interp_PAI_all!1 					= "case"
;var_interp_PAI_all&case 				= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)
;
;; -----------SI index-----------------------------------------------
;var_interp_SI_all!0 					= "case"
;var_interp_SI_all&case 					= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)
;
;; -----------skill score defined by Taylor (2001)-----------------------------------------------
;var_interp_Ra_all!0				 	= "region"
;var_interp_Ra_all!1 					= "case"
;var_interp_Ra_all&case 					= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)
;
;var_interp_ratio_all!0				 	= "region"
;var_interp_ratio_all!1 					= "case"
;var_interp_ratio_all&case 				= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)
;
;var_interp_SS_all!0				 	= "region"
;var_interp_SS_all!1 					= "case"
;var_interp_SS_all&case 					= (/"obs--obs","amip-def","cmip-def","amip-new","cmip-new"/)
;

;do ii=0,dimsizes(latSa)-1

;print("--------------"+regions(ii)+"-----------------")
;; ----------organize all indexes in one table------------------------------------------------
;print("cases		PAI		SI		pat_cor		ratio		skill_score")
;print(var_interp_Ra_all&case+"	"\
;				+sprintf("%9.3f",var_interp_PAI_all(ii,:))+"	"\
;				+sprintf("%9.3f",var_interp_SI_all)+"	"\
;				+sprintf("%9.3f",var_interp_Ra_all(ii,:))+"	"\
;				+sprintf("%9.3f",var_interp_ratio_all(ii,:))+"	"\
;				+sprintf("%9.3f",var_interp_SS_all(ii,:))+"	"\
;     )
;
;end do ; do ii=

end do; do ivar=

                               
end
