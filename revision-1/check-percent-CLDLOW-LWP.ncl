
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

function regrid_diff_time(time[*]:numeric,\
                          var_obs[*][*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                          var_def[*][*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		          var_new[*][*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

lat_spc			 	= fspan(-90,90,72)
lon_spc				= fspan(0,360,144)

;horizontal interpolation
var_interp 			= time_lat_lon_interpolate_3vars_specgrid (dimsizes(time), \
                                			         	   var_obs,lat_obs,lon_obs, \
							 		   var_def,lat_def,lon_def, \
							 	 	   var_new,lat_new,lon_new, \
									   lat_spc,lon_spc)
var_interp_diff 		= var_interp

var_interp_diff 		= var_interp
; 2018-06-18 modify it into percentage value rather than absolute value
var_interp_diff(0,:,:,:) 	= (var_interp(1,:,:,:) - var_interp(0,:,:,:));*100/where(var_interp(0,:,:,:).ne.0,var_interp(0,:,:,:),var_interp@_FillValue) ; def-obs
var_interp_diff(1,:,:,:) 	= (var_interp(2,:,:,:) - var_interp(0,:,:,:));*100/where(var_interp(0,:,:,:).ne.0,var_interp(0,:,:,:),var_interp@_FillValue) ; new-obs
var_interp_diff(2,:,:,:) 	= (var_interp(2,:,:,:) - var_interp(1,:,:,:));*100/where(var_interp(1,:,:,:).ne.0,var_interp(1,:,:,:),var_interp@_FillValue) ; new-obs

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin

CLDLOW_levs			= (/5.,10,15,20,25,30,40,50,60,70,75,80,85,90,95/)
CLDLOW_Diff_levs		= (/-50.,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)

TGCLDLWP_levs			= (/10,25,50,75,100,150,175,200,250,300,350,400,450,500,550/)
TGCLDLWP_Diff_levs		= (/-120,-100,-80,-60,-40,-20,-10,0,10,20,40,60,80,100,120/)

FSNS_levs			= (/0,25,50,75,100,125,150,175,200,220,240,260,280,300,320/)
FSNS_Diff_levs			= (/-100,-75,-50,-40,-30,-20,-10,0,10,20,30,40,50,75,100/)

standard_Diff_levs		= (/-50.,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)


var				=(/"FSNS","CLDLOW","TGCLDLWP"/)
var_levs			= (/FSNS_levs,CLDLOW_levs,TGCLDLWP_levs/)
varDiff_levs			= (/FSNS_Diff_levs,CLDLOW_Diff_levs,TGCLDLWP_Diff_levs/)
;varDiff_levs			= (/standard_Diff_levs,standard_Diff_levs,standard_Diff_levs/)

var_units			= (/"W/m2","fraction","g/m2"/)

compare_type 			= (/"cmip-obs"/)

case				= (/"def","new"/)

Dir				= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
	print(compare_type(icomp))

do icase=0,dimsizes(case)-1
	print(case(icase))

do ivar=0,dimsizes(var)-1
	if(compare_type(icomp) .eq."amip-obs")then
		casename 	= (/"","FC5_f09f09_MG15_COSP","new-cldfrac16-1d-MG15-revise-2-COSP"/)
		fdir 		= (/Dir+"obs/",Dir+"amip/",Dir+"amip/"/)

	end if
	if(compare_type(icomp) .eq."cmip-obs")then
		casename 	= (/"","BC5_f19g16","BC5_f19g16_mac2"/)
		fdir 		= (/Dir+"obs/",Dir+"atm/",Dir+"atm/"/)

	end if

fapend 			= (/"_MONTHS.nc", \
          		    ".cam.h0.MONTHS.nc", \
	  		    ".cam.h0.MONTHS.nc" /)


if(var(ivar).eq."SST")then
	casename(0) 		= "HadISST_PD"
	fdir 			= (/Dir+"obs/", \
	        	            Dir+"ocn/", \
			            Dir+"ocn/"/)
fapend 				= (/"_MONTHS.nc", \
          			    ".pop.h.MONTHS.nc", \
	  			    ".pop.h.MONTHS.nc" /)
end if

if(var(ivar).eq."PRECT")then
	casename(0) 		= "GPCP"
end if

if(var(ivar).eq."CLDLOW".or.var(ivar).eq."CLDMED".or.var(ivar).eq."CLDHGH".or.var(ivar).eq."CLDTOT")then
	casename(0) 		= "CLOUDSAT"
end if

if(var(ivar).eq."FSNS" .or. var(ivar).eq."FLNS" .or. var(ivar).eq."LHFLX" .or. var(ivar).eq."SHFLX" .or. var(ivar).eq."NetSfcHeat")then
	casename(0) 		= "JRA25"
end if

if(var(ivar).eq."TGCLDLWP")then
	casename(0)		= "NVAP"
end if


latS 				= -20
latE 				= -5

lonS				= 230;180
lonE				= 255;230



wks   = gsn_open_wks ("x11", "./figure/check_percent_CLDLOW_"+var(ivar)+"_"+compare_type(icomp)+"_"+latS+"_"+latE+"_"+lonS+"_"+lonE+"_online_revision_1") 

result_obs			= get_var_time(fdir(0),fapend(0),casename(0),var(ivar),"obs")
var_obs				= result_obs[0]
lat_obs				= result_obs[1]
lon_obs				= result_obs[2]
time				= result_obs[3]

result_def			= get_var_time(fdir(1),fapend(1),casename(1),var(ivar),"model")
var_def				= result_def[0]
lat_def				= result_def[1]
lon_def				= result_def[2]


result_new			= get_var_time(fdir(2),fapend(2),casename(2),var(ivar),"model")
var_new				= result_new[0]
lat_new				= result_new[1]
lon_new				= result_new[2]

printVarSummary(var_obs)

result 				= regrid_diff_time(time, \
			 			   var_obs,lat_obs,lon_obs,\
			 			   var_def,lat_def,lon_def,\
			 			   var_new,lat_new,lon_new)

var_interp 			= result[0]
var_interp_diff 		= result[1]

printVarSummary(var_interp_diff)
printMinMax(var_interp_diff,True)


delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])


;;;;; zonal mean

latS 				= -20
latE 				= -5

lonS				= 230;180
lonE				= 255;230


; zonal mean over eastern pacific
var_interp_zm 			= dim_avg_n_Wrap(var_interp(:,:,:,{lonS:lonE}),3)
var_interp_diff_zm 		= dim_avg_n_Wrap(var_interp_diff(:,:,:,{lonS:lonE}),3)

var_interp_zm_rev 		= var_interp_zm(case|:,lat|:,time|:)
var_interp_diff_zm_rev 		= var_interp_diff_zm(case|:,lat|:,time|:)

printVarSummary(var_interp_diff_zm_rev)
printMinMax(var_interp_zm_rev,True)
printMinMax(var_interp_diff_zm_rev,True)

;;;; region average : region C (180-105W, 5S-20S)
lat_here                        = var_interp&lat
rad                             = 4.0*atan(1.0)/180.
clat                            = cos(lat_here*rad)
clat!0                          = "lat"
clat&lat                        = lat_here

var_interp_avg			= wgt_areaave_Wrap(var_interp(:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1,1)
var_interp_diff_avg		= wgt_areaave_Wrap(var_interp_diff(:,:,{latS:latE},{lonS:lonE}),clat({latS:latE}),1,1)

printVarSummary(var_interp_avg)

;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx line plot xxxxxxxxxxxxxxxxxxxxxxxxxxx

plotl				= new(3,graphic)
resl				= True
resl@gsnFrame			= False
resl@gsnDraw			= False
resl@xyLineThicknessF		= 3.0
resl@xyDashPattern		= (/0,0,0/)

resl@tmXMajorGrid                = True
resl@tmXMajorGridThicknessF      = 1.0
resl@tmXMajorGridLineDashPattern = 2

resl@tmYMajorGrid                = True
resl@tmYMajorGridThicknessF      = 1.0 
resl@tmYMajorGridLineDashPattern = 2   

resl@gsnLeftString		= var(ivar)

resl@xyLineColors		= (/"grey","red","blue"/)
plotl(0)			= gsn_csm_y(wks,var_interp_avg(:,:),resl)

resl@xyLineColors		= (/"red","blue","black"/)
plotl(1)			= gsn_csm_y(wks,var_interp_diff_avg(:,:),resl)

;xxxxxxxxxxxx add legends xxxxxxxxxxxxxxxxxxxxx
legend1                         = (/"OBS","DEF","NEW"/)
legend2				= (/"DEF-OBS","NEW-OBS","NEW-DEF"/)

lgres                           = True
;lgres@lgItemType                = resl@xyMarkLineMode
lgres@lgDashIndexes		= resl@xyDashPattern
lgres@lgLineThicknessF 		= resl@xyLineThicknessF

lgres@lgLabelFontHeightF        = 0.08
lgres@vpWidthF                  = 0.15
lgres@vpHeightF                 = 0.1
lgres@lgPerimOn                 = False

lgres@lgLineColors		= (/"grey","red","blue"/)
lbid1                           = gsn_create_legend(wks,dimsizes(legend1),legend1,lgres)
lgres@lgLineColors		= (/"red","blue","black"/)
lbid2                           = gsn_create_legend(wks,dimsizes(legend2),legend2,lgres)

amres                           = True
amres@amJust			= "TopLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = -0.50

annoid1                         = gsn_add_annotation(plotl(0),lbid1,amres)
annoid2                         = gsn_add_annotation(plotl(1),lbid2,amres)

gsn_panel(wks,(/plotl(:1)/),(/1,2/),False)


;*************************
; plotting parameters
;*************************
if(var(ivar).eq."PRECT")then
	cmap1				= read_colormap_file("precip_11lev")
	cmap2                           = read_colormap_file("precip_diff_12lev")
else
	cmap1                           = read_colormap_file("WhiteBlueGreenYellowRed")
	cmap2                           = read_colormap_file("NCV_blu_red")
end if

plot 						= new(6,graphic)

res                  				= True
res@gsnFrame 					= False
res@gsnDraw 					= False

res@cnFillOn         				= True
res@cnLinesOn        				= False

res@lbOrientation    				= "Vertical"
res@pmLabelBarOrthogonalPosF 			= -0.04

res@gsnStringFontHeightF			= 0.025
res@tmXBLabelFontHeightF			= 0.025
res@tmYLLabelFontHeightF			= 0.025

res@gsnCenterStringOrthogonalPosF		= 0.02
res@gsnLeftStringOrthogonalPosF			= 0.02
res@gsnRightStringOrthogonalPosF		= 0.02


res@vpXF             				= 0.12                ; default is 0.2 (aspect ratio)
res@vpYF             				= 0.8                 ; default is 0.8
res@vpHeightF        				= 0.4                 ; default is 0.6
res@vpWidthF         				= 0.75                ; default is 0.6

res@tmXBMode 					= "Explicit"
res@tmXBValues 					= fspan(0,11,12)
res@tmXBLabels 					= (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

res@lbLabelStride    				= 2

if(lonS.gt.180)then
res@gsnCenterString 				= (360.-lonE)+"W-"+(360-lonS)+"W"
else
res@gsnCenterString 				= lonS+"E-"+(360-lonE)+"W"
end if

res@gsnRightString 				= var(ivar)

res@lbTitleOn					= True
res@lbTitleString				= var_units(ivar)


res1						= res
;res1@cnFillPalette        			= cmap1(:170,:)
res1@cnFillPalette        			= cmap1

res1@cnLevelSelectionMode			= "ExplicitLevels"
res1@cnLevels					= var_levs(ivar,:)



if(compare_type(icomp).eq."cmip-amip")then
	string_cli 			= (/"OBS","AMIP_"+case(icase),"Coupled_"+case(icase)/)
else
	string_cli 			= (/"OBS","DEF","NEW"/)
end if


do icli = 0,2
	res1@gsnLeftString 		= string_cli(icli)
	plot(icli) 			= gsn_csm_contour(wks, var_interp_zm_rev(icli,{latS:latE},:), res1) 
end do ; do icli=

;gsn_panel(wks,plot(0:2),(/3,1/),False)
;exit

;------------------------------------------------------------------------------------------------------------
res2					= res
res2@cnFillPalette        		= cmap2


if(compare_type(icomp).eq."cmip-amip")then
	string_diff = (/"AMIP_"+case(icase)+"-OBS",\
	               "Coupled_"+case(icase)+"-OBS",\
		       "Coupled_"+case(icase)+"-AMIP_"+case(icase)/)
else
	string_diff = (/"DEF-OBS","NEW-OBS","NEW-DEF"/)
end if

do idiff = 0,2

 res2@cnLevelSelectionMode		= "ExplicitLevels"
 res2@cnLevels				= varDiff_levs(ivar,:)


	res2@gsnLeftString 		= string_diff(idiff)
	plot(idiff+3) 			= gsn_csm_contour(wks, var_interp_diff_zm_rev(idiff,{latS:latE},:), res2) 

end do ; do idff=

; gsn_panel(wks,plot(3:5),(/3,1/),False)


resP 				= True
;resP@gsnPanelLabelBar		= True
resP@amJust			= "TopLeft"
resP@gsnPanelFigureStrings	= (/"a)","b)","c)","d)"/)

;gsn_panel(wks,(/plot(0),plot(1),plot(3),plot(5)/),(/2,2/),resP)


resP2				= True
resP2@amJust			= "TopLeft"
resP2@gsnPanelFigureStrings	= (/"a)","b)","c)","d)","e)","f)"/)

;gsn_panel(wks,(/plot(0),plot(3),plot(1),plot(4),plot(2),plot(5)/),(/3,2/),resP2)


delete(time)

delete(var_interp)
delete(var_interp_diff)
delete(var_interp_zm)
delete(var_interp_diff_zm)
delete(var_interp_zm_rev)
delete(var_interp_diff_zm_rev)

delete(plot)
delete(fapend)
delete(res)

delete(fdir)

end do ; do ivar=

if(compare_type(icomp) .eq."amip-obs" .or. compare_type(icomp) .eq."cmip-obs")then
	break
end if

end do; do icase=

end do; do icomp=
 
                              
end
