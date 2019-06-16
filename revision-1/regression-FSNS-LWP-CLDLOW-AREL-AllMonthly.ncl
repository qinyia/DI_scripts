
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

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

;var				=(/"FSNS","CLDLOW","TGCLDLWP","AREL"/)
var				= (/"AREL","CLDLOW","TGCLDLWP","FSNS"/)
var_units			= (/"micron","fraction","g/m2","W/m2"/)

compare_type 			= (/"cmip-obs"/)

case				= (/"def","new"/)

;Dir				= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"
Dir				= "/home/lyl/WORK1/lyl_cesm1_2_1/archive/"

do icomp = 0,dimsizes(compare_type)-1
	print(compare_type(icomp))

do icase=0,dimsizes(case)-1
	print(case(icase))


do ivar=0,dimsizes(var)-1
	if(compare_type(icomp) .eq."amip-obs")then
		casename 	= (/"FC5_f09f09_MG15_COSP","new-cldfrac16-1d-MG15-revise-2-COSP"/)
		fdir 		= (/Dir,Dir/)

	end if
	if(compare_type(icomp) .eq."cmip-obs")then
		casename 	= (/"BC5_f19g16","BC5_f19g16_mac2"/)
		fdir 		= (/Dir,Dir/)

	end if

fapend 			= (/\
          		    ".cam.h0.MONTHS.nc", \
	  		    ".cam.h0.MONTHS.nc" /)

fdir_append 			= (/\
          		    "/atm/hist/", \
	  		    "/atm/hist/" /)



latS 				= -20
latE 				= -5

lonS				= 230;180
lonE				= 255;230

;get_var_atm_TimeLevLatLon_addfiles(Dir:string,fdir_append:string,case:string,var:string,type:string)
; fils                    = systemfunc ("ls "+Dir+case+fdir_append+"*.nc") ; file paths

if(var(ivar).ne."AREL")then
	result_def			= get_var_atm_TimeLatLon_addfiles(fdir(0),fdir_append(0),casename(0),var(ivar),"model")
	var_def				= result_def[0]
	lat_def				= result_def[1]
	lon_def				= result_def[2]
	time_def			= result_def[3]
	
	result_new			= get_var_atm_TimeLatLon_addfiles(fdir(1),fdir_append(1),casename(1),var(ivar),"model")
	var_new				= result_new[0]
	lat_new				= result_new[1]
	lon_new				= result_new[2]
else
	result_def			= get_var_atm_TimeLevLatLon_addfiles(fdir(0),fdir_append(0),casename(0),var(ivar),"model")
	var_def_tmp			= result_def[0]
	lat_def				= result_def[1]
	lon_def				= result_def[2]
	time				= result_def[3]
	lev_def				= result_def[4]
	
	result_new			= get_var_atm_TimeLevLatLon_addfiles(fdir(1),fdir_append(1),casename(1),var(ivar),"model")
	var_new_tmp			= result_new[0]
	lat_new				= result_new[1]
	lon_new				= result_new[2]

	; get the column-averaged liquid effective radius.
	; 1. let all zero value to missing value to avoid their weights to get the mean.
	var_def_tmp			= where(var_def_tmp.ne.0,var_def_tmp,var_def_tmp@_FillValue)
	var_new_tmp			= where(var_new_tmp.ne.0,var_new_tmp,var_new_tmp@_FillValue)
	printVarSummary(var_def_tmp)
	; 2. average below 680hPa to make sure they are related to the low clouds.
	var_def				= dim_avg_n_Wrap(var_def_tmp(:,{680:},:,:),1)
	var_new				= dim_avg_n_Wrap(var_new_tmp(:,{680:},:,:),1)

	printVarSummary(var_def)
	printMinMax(var_def,True)
	printMinMax(var_new,True)
end if

printVarSummary(var_def)
printVarSummary(var_new)

result 				= regrid_diff_time(time, \
			 			   var_def,lat_def,lon_def,\
			 			   var_new,lat_new,lon_new)

var_interp 			= result[0]
var_interp_diff 		= result[1]

dsizes				= dimsizes(var_interp)
dsizes_diff			= dimsizes(var_interp_diff)

if(ivar.eq.0)then
	varall				= new((/dimsizes(var),dsizes(0),dsizes(1),dsizes(2),dsizes(3)/),float)
	varall!0			= "vars"
	varall&vars			= fspan(0,dimsizes(var)-1,dimsizes(var))

	varall_diff			= new((/dimsizes(var),dsizes_diff(0),dsizes_diff(1),dsizes_diff(2)/),float)
	varall_diff!0			= "vars"
	varall_diff&vars		= fspan(0,dimsizes(var)-1,dimsizes(var))

	varall(ivar,:,:,:,:)		= var_interp
	varall_diff(ivar,:,:,:)		= var_interp_diff
else
	varall(ivar,:,:,:,:)		= var_interp
	varall_diff(ivar,:,:,:)		= var_interp_diff
end if

delete([/var_def,lat_def,lon_def,var_new,lat_new,lon_new/])
delete([/dsizes,dsizes_diff/])
;delete([/var_def_tmp,var_new_tmp/])

end do ; do ivar=

printVarSummary(varall)

;;;;; zonal mean

latS 				= -20
latE 				= -5

lonS				= 230;180
lonE				= 255;230


; zonal mean over eastern pacific
varall_zm 			= dim_avg_n_Wrap(varall(:,:,:,:,{lonS:lonE}),4)
varall_diff_zm 			= dim_avg_n_Wrap(varall_diff(:,:,:,{lonS:lonE}),3)

varall_zm_rev 			= varall_zm(vars|:,case|:,lat|:,time|:)
varall_diff_zm_rev 		= varall_diff_zm(vars|:,lat|:,time|:)

;;;; region average : region C (180-105W, 5S-20S)
lat_here                        = varall&lat
rad                             = 4.0*atan(1.0)/180.
clat                            = cos(lat_here*rad)
clat!0                          = "lat"
clat&lat                        = lat_here

varall_avg			= wgt_areaave_Wrap(varall(:,:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1,1)
varall_diff_avg			= wgt_areaave_Wrap(varall_diff(:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1,1)

printVarSummary(varall_avg)

AREL_avg			= varall_avg(0,:,:)
CLDLOW_avg			= varall_avg(1,:,:)
LWP_avg				= varall_avg(2,:,:)

; in-cloud LWP
incld_LWP_avg			= varall_avg(2,:,:)
incld_LWP_avg			= varall_avg(2,:,:)*100./where(varall_avg(1,:,:).ne.0,varall_avg(1,:,:),varall_avg@_FillValue)
incld_LWP_avg_diff		= incld_LWP_avg(0,:)
incld_LWP_avg_diff		= incld_LWP_avg(1,:)-incld_LWP_avg(0,:)
printVarSummary(incld_LWP_avg_diff)
printMinMax(incld_LWP_avg_diff,True)

; optical depth \approx = LWP/Re
; from micron to meter
tau_avg				= varall_avg(2,:,:)
tau_avg				= varall_avg(2,:,:)/where(varall_avg(0,:,:).ne.0,varall_avg(0,:,:),varall_avg@_FillValue)
tau_avg_diff			= tau_avg(0,:)
tau_avg_diff			= tau_avg(1,:)-tau_avg(0,:)


dsizes				= dimsizes(varall_avg)
dsizes_diff			= dimsizes(varall_diff_avg)
data				= new((/dimsizes(var)+2,dsizes(1),dsizes(2)/),float)
data(:dimsizes(var)-1,:,:)	= varall_avg
data(dimsizes(var),:,:)		= incld_LWP_avg
data(dimsizes(var)+1,:,:)	= tau_avg

data_diff				= new((/dimsizes(var)+2,dsizes_diff(1)/),float)
data_diff(:dimsizes(var)-1,:)		= varall_diff_avg
data_diff(dimsizes(var),:)		= incld_LWP_avg_diff
data_diff(dimsizes(var)+1,:)		= tau_avg_diff

printVarSummary(data_diff)
printMinMax(data_diff,True)

delete(dsizes)
delete(dsizes_diff)

;=====================================================================================================
; Qinyi
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; regression 
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

; normalized by def abolute value: get relative change percentage
data_diff_percent			= data_diff/where(data(:,0,:).ne.0,data(:,0,:),data_diff@_FillValue)
copy_VarCoords(data_diff,data_diff_percent)
data_pro				= data_diff_percent
printVarSummary(data_pro)

; get anomaly
data_diff_ano				= data_pro
data_diff_reg				= data_diff_ano
dsizes					= dimsizes(data_pro)
maxx_all				= new((/dsizes(0)/),float)
avg_all					= new((/dsizes(0)/),float)

;=====================================================================================================
;xxxxxxxxxxxxxxxxxxxxxxx Parameters Settings xxxxxxxxxxxxxxx
newlist					= (/"AREL","CLDLOW","TGCLDLWP","FSNS","incld-LWP","tau"/)
					;      0	1	2	3	4	5

typestring				= "raw";"normalized"
w_intercept				= False

plottype				= "x11"
timetype				= "AllMonthly"
output1					= "./figure/regression-plot-"+timetype+"-"+typestring+"-interpcet-"+w_intercept
output2					= "./figure/regression-FSNS-LWP-CLDLOW-AREL"+timetype+"-"+typestring+"-interpcet-"+w_intercept
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

do ivar=0,dsizes(0)-1

avg_all(ivar)					= dim_avg_Wrap(data_pro(ivar,:))
data_diff_ano(ivar,:)				= data_pro(ivar,:)-avg_all(ivar)
maxx_all(ivar)					= max(abs(data_diff_ano(ivar,:)))

if(typestring.eq."normalized")then
	; normalized by "max"
	data_diff_reg(ivar,:)			= data_diff_ano(ivar,:)/where(maxx_all(ivar).ne.0,maxx_all(ivar),data@_FillValue)
else if(typestring.eq."raw")then
	; no normalized
	data_diff_reg(ivar,:)			= data_pro(ivar,:)
end if
end if

end do ; do ivar=

printVarSummary(data_diff_reg)
printMinMax(data_diff_reg,True)

regy					= data_diff_reg(3,:)
regx1_tmp				= data_diff_reg((/0/),:) ; AREL
regx1					= regx1_tmp

regx2_tmp				= data_diff_reg((/1/),:) ; CLDLOW
regx2					= regx2_tmp

regx3_tmp				= data_diff_reg((/4/),:) ; incldLWP
regx3					= regx3_tmp

regx4_tmp				= data_diff_reg((/2/),:) ; gridmean LWP
regx4					= regx4_tmp

regx5_tmp				= data_diff_reg((/0,1,4/),:) ; AREL,CLDLOW,incldLWP
regx5					= regx5_tmp(time|:,vars|:)


opt 					= True
opt@print_anova 			= False
opt@print_data  			= False
b1   					=  reg_multlin_stats(regy,regx1,opt)
b2   					=  reg_multlin_stats(regy,regx2,opt)
b3   					=  reg_multlin_stats(regy,regx3,opt)
b4  					=  reg_multlin_stats(regy,regx4,opt)
b5  					=  reg_multlin_stats(regy,regx5,opt)

yest1					= b1(0)+regx1(:)*b1(1)
yest2					= b2(0)+regx2(:)*b2(1)
yest3					= b3(0)+regx3(:)*b3(1)
yest4					= b4(0)+regx4(:)*b4(1)
yest5					= b5(0)+regx5(:,0)*b5(1)+regx5(:,1)*b5(2)+regx5(:,2)*b5(3)

;xxxxxxxxxxxxxxxxxxxxxxxxxxxx
; this is for testing consistency between normalized and non-normalized.--------Yes, this works well. 
if(typestring.eq."normalized")then
	print("yest5="+sprintf("%6.3f",b5(0)-maxx_all(3)*(b5(1)*avg_all(0)/maxx_all(0)+b5(2)*avg_all(1)/maxx_all(1)+b5(3)*avg_all(4)/maxx_all(4))+avg_all(3))+"+"+sprintf("%6.3f",b5(1)*maxx_all(3)/maxx_all(0))+"*AREL+"+sprintf("%6.3f",b5(2)*maxx_all(3)/maxx_all(1))+"*CLDLOW+"+sprintf("%6.3f",b5(3)*maxx_all(3)/maxx_all(4))+"*incld-LWP")
else
	print("yest5="+sprintf("%6.3f",b5(0))+"+"+sprintf("%6.3f",b5(1))+"*AREL+"+sprintf("%6.3f",b5(2))+"*CLDLOW+"+sprintf("%6.3f",b5(3))+"*incld-LWP")
end if

;xxxxxxxxxxxxxxxxxxxxxxxxxxxx
; estimate the contribution from incld-LWP and AREL to Total cloud fraction percentage change.
if(typestring.eq."raw")then
if(w_intercept)then
	;					           AREL             incld-LWP         CF
	total_increase_percent_CF		= (b5(1)*avg_all(0)+b5(3)*avg_all(4)+b5(0))/b5(2)+avg_all(1)
else
	total_increase_percent_CF		= (b5(1)*avg_all(0)+b5(3)*avg_all(4))/b5(2)+avg_all(1)
	;print(total_increase_percent_CF)
end if
end if
;xxxxxxxxxxxxxxxxxxxxxxxxxxxx

data_reg				= (/regy,yest1,yest2,yest3,yest4,yest5/)

;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx regression plot xxxxxxxxxxxxxxxxxxxxxxxxxxx
wks					= gsn_open_wks(plottype,output1)
resg					= True

resg@gsnFrame				= False
resg@gsnDraw				= False

resg@vpWidthF				= 0.7
resg@vpHeightF				= 0.4

resg@xyMarkLineModes			= (/"MarkLines","Lines","Lines","Lines","Lines","Lines"/)

resg@xyMarkers				= (/16,-1,-1,-1,-1,-1/)
resg@xyMarkerColors			= (/"red","white","white","white","white","white"/)
resg@xyMarkerSizeF			= 0.01
resg@xyDashPatterns			= (/0,0,0,0,0,0/)
resg@xyLineColors			= (/"red","blue","green","orange","purple","black"/)
resg@xyLineThicknessF			= 3.0

resg@tiMainString			= "Relative to "+typestring+" simulation"

resg@gsnRightString			= "total_increase_percent_CF= "+total_increase_percent_CF

plotg					= gsn_csm_y(wks,data_reg,resg)

;xxxxxxxxxxxx add legends xxxxxxxxxxxxxxxxxxxxx
legend1                         = (/"regy",\
				    "yest1="+sprintf("%6.3f",b1(0))+"+"+sprintf("%6.3f",b1(1))+"*AREL",\
				    "yest2="+sprintf("%6.3f",b2(0))+"+"+sprintf("%6.3f",b2(1))+"*CLDLOW",\
				    "yest3="+sprintf("%6.3f",b3(0))+"+"+sprintf("%6.3f",b3(1))+"*incld-LWP",\
				    "yest4="+sprintf("%6.3f",b4(0))+"+"+sprintf("%6.3f",b4(1))+"*LWP",\
				    "yest5="+sprintf("%6.3f",b5(0))+"+"+sprintf("%6.3f",b5(1))+"*AREL+"+sprintf("%6.3f",b5(2))+"*CLDLOW+"+sprintf("%6.3f",b5(3))+"*incld-LWP"\
				  /)

lgres                           = True
lgres@lgMonoItemType		= False
lgres@lgItemTypes               = resg@xyMarkLineModes
lgres@lgDashIndexes		= resg@xyDashPatterns
lgres@lgMarkerIndexes		= resg@xyMarkers
lgres@lgMarkerColors		= resg@xyMarkerColors
lgres@lgLineThicknessF	 	= resg@xyLineThicknessF
lgres@lgLineColors		= resg@xyLineColors

lgres@lgLabelFontHeightF        = 0.2
lgres@vpWidthF                  = 0.7
lgres@vpHeightF                 = 0.3
lgres@lgPerimOn                 = False

lbid1                           = gsn_create_legend(wks,dimsizes(legend1),legend1,lgres)

amres                           = True
amres@amJust			= "TopLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = 0.70

annoid1                         = gsn_add_annotation(plotg,lbid1,amres)

draw(plotg)
frame(wks)

delete(legend1)
delete(amres)
delete(lgres)
delete(annoid1)
delete(lbid1)

exit

;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx line plot xxxxxxxxxxxxxxxxxxxxxxxxxxx


wks				= gsn_open_wks(plottype,output2)
plotl				= new((/3,dimsizes(newlist)/),graphic)

do ivar=0,dimsizes(newlist)-1

resl				= True
resl@gsnFrame			= False
resl@gsnDraw			= False

resl@vpWidthF			= 0.7
resl@vpHeightF			= 0.4

resl@xyLineThicknessF		= 3.0
resl@xyDashPattern		= (/0,0,0/)

;resl@tmXMajorGrid                = True
;resl@tmXMajorGridThicknessF      = 1.0
;resl@tmXMajorGridLineDashPattern = 2
;
;resl@tmYMajorGrid                = True
;resl@tmYMajorGridThicknessF      = 1.0 
;resl@tmYMajorGridLineDashPattern = 2   

resl@tmXBLabelFontHeightF        = 0.04
resl@tmYLLabelFontHeightF        = 0.04
resl@gsnStringFontHeightF	 = 0.04

resl@gsnLeftString		= newlist(ivar)

resl@trYMaxF			= max(data(ivar,:,:))*1.2

resl@xyLineColors		= (/"red","blue"/)
plotl(0,ivar)			= gsn_csm_y(wks,data(ivar,:,:),resl)

resl@xyLineColors		= (/"grey","red"/)
; absolute difference value
plotl(1,ivar)			= gsn_csm_y(wks,data_diff(ivar,:),resl)

; percentage
plotl(2,ivar)			= gsn_csm_y(wks,data_diff_reg(ivar,:),resl)

;xxxxxxxxxxxx add legends xxxxxxxxxxxxxxxxxxxxx
legend1                         = (/"DEF","NEW"/)
legend2				= (/"NEW-DEF"/)

lgres                           = True
;lgres@lgItemType                = resl@xyMarkLineMode
lgres@lgDashIndexes		= resl@xyDashPattern
lgres@lgLineThicknessF 		= resl@xyLineThicknessF

lgres@lgLabelFontHeightF        = 0.15
lgres@vpWidthF                  = 0.3
lgres@vpHeightF                 = 0.2
lgres@lgPerimOn                 = False

lgres@lgLineColors		= (/"red","blue"/)
lbid1                           = gsn_create_legend(wks,dimsizes(legend1),legend1,lgres)
lgres@lgLineColors		= (/"grey","red"/)
lbid2                           = gsn_create_legend(wks,dimsizes(legend2),legend2,lgres)

amres                           = True
amres@amJust			= "TopLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = -0.50

annoid1                         = gsn_add_annotation(plotl(0,ivar),lbid1,amres)
annoid2                         = gsn_add_annotation(plotl(1,ivar),lbid2,amres)

;gsn_panel(wks,(/plotl(:1,ivar)/),(/1,2/),False)


delete([/annoid1,annoid2,lbid1,lbid2/])

end do ; do ivar=

resP				= True
resP@gsnPanelRight      	= 0.95
resP@gsnPanelLeft		= 0.05
resP@gsnPanelYWhiteSpacePercent = 5.0
resP@gsnPanelXWhiteSpacePercent = 5.0

resP@gsnPanelMainString		= "Relative to "+typestring+" simulation"


plot1				= ndtooned((/plotl(0,:3),plotl(1,:3),plotl(2,:3)/))
gsn_panel(wks,plot1,(/3,dimsizes(newlist(:3))/),resP)

plot2				= ndtooned((/plotl(0,4:),plotl(1,4:),plotl(2,4:)/))
gsn_panel(wks,plot2,(/3,dimsizes(newlist(4:))/),resP)


if(compare_type(icomp) .eq."amip-obs" .or. compare_type(icomp) .eq."cmip-obs")then
	break
end if

end do; do icase=

end do; do icomp=
 
                              
end
