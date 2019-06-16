load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "/home/lyl/WORK3/qinyi/double-ITCZ/functions_lat_interpolate.ncl"
load "/home/lyl/WORK3/qinyi/double-ITCZ/functions_get_resource.ncl"
load "/home/lyl/WORK3/qinyi/double-ITCZ/functions_get_var.ncl"

load "./taylor_stats.ncl"

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ttest_2D(xtmp[*][*]:numeric,ytmp[*][*]:numeric)
begin

xAve = dim_avg (xtmp)              ; calculate means at each grid point 
yAve = dim_avg (ytmp)
xVar = dim_variance (xtmp)         ; calculate variances
yVar = dim_variance (ytmp)
print(xAve+"		"+xVar+"	"+yAve+"	"+yVar)
exit

sigr = 0.05                        ; critical sig lvl for r

xN   = dimsizes(xtmp(0,:))
yN   = dimsizes(ytmp(0,:))

xN   	= 12
yN	= 12 
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

;print("qinyi "+xVar+"	"+yVar)
;exit

sigr = 0.05                        ; critical sig lvl for r

xN   = dimsizes(xtmp)
yN   = dimsizes(ytmp)

xN   	= 12;50;80
yN	= 12;50;80

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

fdir = (/"/home/share3/lyl/work3/qinyi/double-ITCZ/data/obs/", \
         "/home/share3/lyl/work3/qinyi/double-ITCZ/data/amip/", \
	 "/home/share3/lyl/work3/qinyi/double-ITCZ/data/atm/"/)

do ivar=0,dimsizes(var)-1


plot = new((/dimsizes(season),6/),graphic)


fapend 			= (/"_MONTHS.nc", \
          		    ".cam.h0.MONTHS.nc", \
	  		    ".cam.h0.MONTHS.nc" /)



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


result_obs 			= get_var_time(fdir(0),fapend(0),casename(0),var(ivar),"obs")
var_obs 			= result_obs[0]
lat_obs   			= result_obs[1]
lon_obs   			= result_obs[2]

result_amip 			= get_var_time(fdir(1),fapend(1),casename(1),var(ivar),"model")
var_amip 			= result_amip[0]
lat_amip   			= result_amip[1]
lon_amip   			= result_amip[2]
time				= result_amip[3]

result_cmip 			= get_var_time(fdir(2),fapend(2),casename(2),var(ivar),"model")
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

var_interp_PAI_all!0		 	= "region"
var_interp_PAI_all&region		= fspan(0,dimsizes(lonSa)-1,dimsizes(lonSa))
var_interp_PAI_all!1			= "case"
var_interp_PAI_all&case			= fspan(0,4,5)
var_interp_PAI_all!2			= "time"
var_interp_PAI_all&time			= fspan(0,ntime-1,ntime)


printVarSummary(var_interp_PAI_all)
;print(var_interp_PAI_all(0,2,:))

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
var_interp_S_all				= (1+var_interp_Ra_all)^4/(4.*(var_interp_ratio_all+1./var_interp_ratio_all)^2)

var_interp_S_all!0		 	= "region"
var_interp_S_all&region			= fspan(0,dimsizes(lonSa)-1,dimsizes(lonSa))
var_interp_S_all!1			= "case"
var_interp_S_all&case			= fspan(0,4,5)
var_interp_S_all!2			= "time"
var_interp_S_all&time			= fspan(0,ntime-1,ntime)

end do; do icase=
printVarSummary(var_interp_S_all)
print(var_interp_S_all(0,1,:))

; ttest to identify whether two simulations show very significant difference.



x_PAI			= var_interp_PAI_all(:,2,:)-var_interp_PAI_all(:,0,:)
y_PAI			= var_interp_PAI_all(:,4,:)-var_interp_PAI_all(:,0,:)
results			= ttest_2D(x_PAI,y_PAI)
prob_PAI		= results[0]
probf_PAI		= results[1]
print(dim_avg(x_PAI)+"	"+dim_avg(y_PAI)+"		"+prob_PAI+"	"+probf_PAI)
exit

x_PAI			= var_interp_PAI_all(1,2,:)-var_interp_PAI_all(1,0,:)
y_PAI			= var_interp_PAI_all(1,4,:)-var_interp_PAI_all(1,0,:)
results			= ttest_1D(x_PAI,y_PAI)
prob_PAI		= results[0]
probf_PAI		= results[1]
print(dim_avg(x_PAI)+"	"+dim_avg(y_PAI)+"		"+prob_PAI+"	"+probf_PAI)
exit
exit



x_S			= var_interp_S_all(:,2,:)-var_interp_S_all(:,0,:)
y_S			= var_interp_S_all(:,4,:)-var_interp_S_all(:,0,:)
results			= ttest_2D(x_S,y_S)
prob_S			= results[0]
probf_S			= results[1]
print(dim_avg(x_S)+"	"+dim_avg(y_S)+"		"+prob_S+"	"+probf_S)
exit




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
exit


end do; do ivar=

                               
end
