load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "/home/lyl/WORK3/qinyi/double-ITCZ/functions_lat_interpolate.ncl"
load "/home/lyl/WORK3/qinyi/double-ITCZ/functions_get_resource.ncl"
load "/home/lyl/WORK3/qinyi/double-ITCZ/functions_get_var.ncl"

load "./taylor_stats.ncl"

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff(var_obs[*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                     var_amip[*][*]:numeric,lat_amip[*]:numeric,lon_amip[*]:numeric,\
		     var_cmip[*][*]:numeric,lat_cmip[*]:numeric,lon_cmip[*]:numeric)
begin
;horizontal interpolation
var_interp = lat_lon_interpolate_3vars(var_obs,lat_obs,lon_obs, \
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


result_obs 			= get_var_2D(fdir(0),fapend(0),casename(0),var(ivar),"obs")
var_obs 			= result_obs[0]
lat_obs   			= result_obs[1]
lon_obs   			= result_obs[2]

result_amip 			= get_var_2D(fdir(1),fapend(1),casename(1),var(ivar),"model")
var_amip 			= result_amip[0]
lat_amip   			= result_amip[1]
lon_amip   			= result_amip[2]

result_cmip 			= get_var_2D(fdir(2),fapend(2),casename(2),var(ivar),"model")
var_cmip 			= result_cmip[0]
lat_cmip   			= result_cmip[1]
lon_cmip   			= result_cmip[2]


result 				= regrid_diff(var_obs,lat_obs,lon_obs,\
	             			      var_amip,lat_amip,lon_amip,\
		     			      var_cmip,lat_cmip,lon_cmip)
var_interp 			= result[0] ; obs,amip,cmip
var_interp_diff 		= result[1] ; obs,amip,cmip

printVarSummary(var_interp)			; (case,lat,lon) = (3, 72, 144)

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
	var_interp_NH_all 				= new((/dimsizes(lonSa),5/),float)
	var_interp_SH_all 				= new((/dimsizes(lonSa),5/),float)
	var_interp_tropical_all 			= new((/dimsizes(lonSa),5/),float)
	var_interp_PAI_all 				= new((/dimsizes(lonSa),5/),float)
end if

do ii=0,dimsizes(lonSa)-1
var_interp_NH 			= wgt_areaave_Wrap(var_interp(:,{0:20},{lonSa(ii):lonEa(ii)}),  clat({0:20}),1.0,1)
var_interp_SH 			= wgt_areaave_Wrap(var_interp(:,{-20:0},{lonSa(ii):lonEa(ii)}), clat({-20:0}),1.0,1)
var_interp_tropical 		= wgt_areaave_Wrap(var_interp(:,{-20:20},{lonSa(ii):lonEa(ii)}),clat({-20:20}),1.0,1)
var_interp_PAI 			= (var_interp_NH-var_interp_SH)/var_interp_tropical

;;; this varible is used to put all results from two cases.
if(icase.eq.0)then
	var_interp_NH_all(ii,3*icase:3*icase+2) 		= var_interp_NH
	var_interp_SH_all(ii,3*icase:3*icase+2) 		= var_interp_SH
	var_interp_tropical_all(ii,3*icase:3*icase+2) 		= var_interp_tropical
	var_interp_PAI_all(ii,3*icase:3*icase+2) 		= var_interp_PAI
else
	var_interp_NH_all(ii,3*icase:3*icase+1) 		= var_interp_NH(1:2)
	var_interp_SH_all(ii,3*icase:3*icase+1) 		= var_interp_SH(1:2)
	var_interp_tropical_all(ii,3*icase:3*icase+1) 		= var_interp_tropical(1:2)
	var_interp_PAI_all(ii,3*icase:3*icase+1) 		= var_interp_PAI(1:2)
end if
delete(var_interp_NH)
delete(var_interp_SH)
delete(var_interp_tropical)
delete(var_interp_PAI)

end do ; do ii=

;-------------------------------------------
; calculate southern ITCZ (SI) index: 
; annual mean precipitation over (20S-0, 100-150W)
;-------------------------------------------
latS				= -20
latE				= 0
lonS				= 210
lonE				= 260
var_interp_SI 			= wgt_areaave(var_interp(:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1.0,1)

;;; this varible is used to put all results from two cases.
if(icase.eq.0)then
	var_interp_SI_all 				= new((/5/),float)
	var_interp_SI_all(3*icase:3*icase+2) 		= var_interp_SI
else
	var_interp_SI_all(3*icase:3*icase+1) 		= var_interp_SI(1:2)
end if
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
var_interp_Ra_all 		= new((/dimsizes(latSa),5/),float)
var_interp_ratio_all		= new((/dimsizes(latSa),5/),float)
end if


do ii=0,dimsizes(latSa)-1

pr_obs				= var_interp(0,{latSa(ii):latEa(ii)},{lonSa(ii):lonEa(ii)})
pr_amip				= var_interp(1,{latSa(ii):latEa(ii)},{lonSa(ii):lonEa(ii)})
pr_cmip				= var_interp(2,{latSa(ii):latEa(ii)},{lonSa(ii):lonEa(ii)})

; try to directly use the function from Taylor diagram: taylor_stats
tstVar				= pr_obs
refVar				= pr_obs
stat_taylor			= taylor_stats(tstVar,refVar,clat({latSa(ii):latEa(ii)}),0)
Ra_obs				= stat_taylor(0) ;correlation
ratio_obs			= stat_taylor(1) ; ratio

tstVar				= pr_amip
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
printVarSummary(var_interp_S_all)


; -------------get Root Mean square error/deviation (RMSE)
if(icase.eq.0)then
var_interp_rmse_all		= new((/5/),float)
end if

wgt				= conform(var_interp(0,:,:),clat,0)
copy_VarCoords_1(var_interp(lat|:,lon|:,case|:),wgt)
printVarSummary(wgt)

lat1				= -30.
lat2				= 30.
lon1				= 0.
lon2				= 360.

if(icase.eq.0)then
var_interp_rmse_all(3*icase)		= wgt_arearmse2(var_interp(0,{lat1:lat2},{lon1:lon2}),var_interp(0,{lat1:lat2},{lon1:lon2}),wgt({lat1:lat2},{lon1:lon2}),0)
var_interp_rmse_all(3*icase+1)		= wgt_arearmse2(var_interp(0,{lat1:lat2},{lon1:lon2}),var_interp(1,{lat1:lat2},{lon1:lon2}),wgt({lat1:lat2},{lon1:lon2}),0)
var_interp_rmse_all(3*icase+2)		= wgt_arearmse2(var_interp(0,{lat1:lat2},{lon1:lon2}),var_interp(2,{lat1:lat2},{lon1:lon2}),wgt({lat1:lat2},{lon1:lon2}),0)
else
var_interp_rmse_all(3*icase)		= wgt_arearmse2(var_interp(0,{lat1:lat2},{lon1:lon2}),var_interp(1,{lat1:lat2},{lon1:lon2}),wgt({lat1:lat2},{lon1:lon2}),0)
var_interp_rmse_all(3*icase+1)		= wgt_arearmse2(var_interp(0,{lat1:lat2},{lon1:lon2}),var_interp(2,{lat1:lat2},{lon1:lon2}),wgt({lat1:lat2},{lon1:lon2}),0)
end if

end do; do icase=
print(var_interp_rmse_all)
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
;=============================Make a table=================================
; Header
  ncr1  = (/1,1/)               ; 1 row, 1 column
;  x1    = (/0.005,0.335/)       ; Start and end X

  x1    = (/0.005,0.100/)       ; Start and end X
  y1    = (/0.900,0.995/)       ; Start and end Y
  text1 = "CASES"

; Sub header
  ncr2  = (/1,17/)               ; 1 rows, 18 columns
  x2    = (/0.100,0.995/)       ; Start and end X
  y2    = (/0.900,0.995/)       ; Start and end Y
  text2 = (/"PAI_GL","PAI_PA","PAI_SA","PAI_AI",\
  	    "SI_PA",\
	    "Ra_GL","Ra_PA","Ra_SA","Ra_AI",\
	    "SD_GL","SD_PA","SD_SA","SD_AI",\
	    "S_GL","S_PA","S_SA","S_AI"/)

; cases index
ncr3	= (/5,1/)
x3	= (/0.005,0.100/)
y3	= (/0.005,0.900/)

text3_tmp	= (/"obs","AMIP-def","CMIP-def","AMIP-new","CMIP-new"/)
text3		= transpose(text3_tmp)

; Main table body
  ncr4  = (/5,17/)              ; 5 rows, 17 columns
  x4    = (/0.100,0.995/)       ; Start and end X
  y4    = (/0.005,0.900/)       ; Start and end Y

  text4_tmp		= (/ \
				sprintf("%9.3f",var_interp_PAI_all(0,:)),\
				sprintf("%9.3f",var_interp_PAI_all(1,:)),\
				sprintf("%9.3f",var_interp_PAI_all(2,:)),\
				sprintf("%9.3f",var_interp_PAI_all(3,:)),\
				sprintf("%9.3f",var_interp_SI_all),\
				sprintf("%9.3f",var_interp_Ra_all(0,:)),\
				sprintf("%9.3f",var_interp_Ra_all(1,:)),\
				sprintf("%9.3f",var_interp_Ra_all(2,:)),\
				sprintf("%9.3f",var_interp_Ra_all(3,:)),\
				sprintf("%9.3f",var_interp_ratio_all(0,:)),\
				sprintf("%9.3f",var_interp_ratio_all(1,:)),\
				sprintf("%9.3f",var_interp_ratio_all(2,:)),\
				sprintf("%9.3f",var_interp_ratio_all(3,:)),\
				sprintf("%9.3f",var_interp_S_all(0,:)),\
				sprintf("%9.3f",var_interp_S_all(1,:)),\
				sprintf("%9.3f",var_interp_S_all(2,:)),\
				sprintf("%9.3f",var_interp_S_all(3,:))\
		  /)

text4				= transpose(text4_tmp)


tmp_data_tmp		= (/ \
				var_interp_PAI_all(0,:),\
				var_interp_PAI_all(1,:),\
				var_interp_PAI_all(2,:),\
				var_interp_PAI_all(3,:),\
				var_interp_SI_all,\
				var_interp_Ra_all(0,:),\
				var_interp_Ra_all(1,:),\
				var_interp_Ra_all(2,:),\
				var_interp_Ra_all(3,:),\
				var_interp_ratio_all(0,:),\
				var_interp_ratio_all(1,:),\
				var_interp_ratio_all(2,:),\
				var_interp_ratio_all(3,:),\
				var_interp_S_all(0,:),\
				var_interp_S_all(1,:),\
				var_interp_S_all(2,:),\
				var_interp_S_all(3,:)\
		  /)
tmp_data		= transpose(tmp_data_tmp)
printVarSummary(tmp_data) ; (5,17)


tmp_data_diff					= tmp_data
do irow=0,4
tmp_data_diff(irow,:)			= tmp_data(irow,:)-tmp_data(0,:)
end do

;text4			= sprintf("%9.3f",tmp_data_diff)

colors = (/"Red","Green","White"/)

; Assign a color to a particular range of data.
  cols_for_vals = new(dimsizes(tmp_data),"string")

do icol=0,16
	maxval_a		= max((/abs(tmp_data_diff(1,icol)),abs(tmp_data_diff(3,icol))/))
do irow=1,3,2 ;amip
	cols_for_vals(irow,icol)= where(abs(tmp_data_diff(irow,icol)).lt.maxval_a,"Green","Red")	
end do ; do irow=

	maxval_c		= max((/abs(tmp_data_diff(2,icol)),abs(tmp_data_diff(4,icol))/))
do irow=2,4,2 ;cmip
	cols_for_vals(irow,icol)= where(abs(tmp_data_diff(irow,icol)).lt.maxval_c,"Green","Red")	
end do ; do irow=
	
	cols_for_vals(0,icol)	= "white"

delete(maxval_a)
delete(maxval_c)
end do ; do icol=

  wks = gsn_open_wks("x11","./figure/All-DI-index-table")  ; send graphics to PNG file

;
; Main header.
;
  res1               = True
  res1@txFontHeightF = 0.01
  res1@gsFillColor   = "CornFlowerBlue"

  gsn_table(wks,ncr1,x1,y1,text1,res1)

;
; Sub header
;
  res2               = True
  res2@txFontHeightF = 0.01
  res2@gsFillColor   = "Gray"

  gsn_table(wks,ncr2,x2,y2,text2,res2)

  res3 = True             ; Set up resource list
  res3@txFontHeightF = 0.01

  gsn_table(wks,ncr3,x3,y3,text3,res3)

;
; Main body of table.
;
  res4 = True             ; Set up resource list
  res4@txFontHeightF 	= 0.01
  res4@gsFillColor	= cols_for_vals

  gsn_table(wks,ncr4,x4,y4,text4,res4)

  frame(wks)      ; Advance the frame.

exit

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

end do ; do is=


end do; do ivar=

                               
end
