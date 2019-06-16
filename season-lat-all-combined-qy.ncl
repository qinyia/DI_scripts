
load "./functions_lat_interpolate.ncl"
load "./functions_get_var.ncl"
load "./functions_get_resource.ncl"

function regrid_diff_time(time[*]:numeric,\
                          var_obs[*][*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                          var_def[*][*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		          var_new[*][*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin

;horizontal interpolation
var_interp 			= time_lat_lon_interpolate_3vars(dimsizes(time), \
                                			         var_obs,lat_obs,lon_obs, \
							 	var_def,lat_def,lon_def, \
							 	var_new,lat_new,lon_new)
var_interp_diff 		= var_interp

var_interp_diff 		= var_interp
var_interp_diff(0,:,:,:) 	= var_interp(1,:,:,:) - var_interp(0,:,:,:) ; def-obs
var_interp_diff(1,:,:,:) 	= var_interp(2,:,:,:) - var_interp(0,:,:,:) ; new-obs
var_interp_diff(2,:,:,:) 	= var_interp(2,:,:,:) - var_interp(1,:,:,:) ; new-obs

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin

PRECT_levs			= (/.2,.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
CLDLOW_levs			= (/5.,10,15,20,25,30,40,50,60,70,75,80,85,90,95/)
SST_levs			= (/-1,0,1,3,6,9,12,15,18,20,22,24,26,28,29/)

;PRECT_Diff_levs			= (/-8,-6,-4,-3,-2,-1,-.5,0,.5,1,2,3,4,6,8/)
PRECT_Diff_levs			= (/-18,-14,-10,-6,-4,-2,-1,0,1,2,4,6,10,14,18/)
CLDLOW_Diff_levs		= (/-50.,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)
;SST_Diff_levs			= (/-5,-4,-3,-2,-1,-.5,-.2,0,.2,.5,1,2,3,4,5/)
SST_Diff_levs			= (/-6,-5,-4,-3,-2,-1.,-0.5,0,.5,1,2,3,4,5,6/)


;var				=(/"SST","PRECT","CLDTOT","CLDHGH","CLDMED","CLDLOW","FSNS","NetSfcHeat"/)
;var				=(/"NetSfcHeat"/)
;var				=(/"FSNS","FLNS","LHFLX","SHFLX"/)
;var				=(/"FLNS"/)
;var				=(/"LHFLX"/)
;var				=(/"CLDTOT"/)

var				=(/"SST","PRECT"/)
var_levs			= (/SST_levs,PRECT_levs/)
varDiff_levs			= (/SST_Diff_levs,PRECT_Diff_levs/)

;var				=(/"CLDLOW","SST"/)
;var_levs			= (/CLDLOW_levs,SST_levs/)
;varDiff_levs			= (/CLDLOW_Diff_levs,SST_Diff_levs/)

;compare_type 			= (/"amip-obs","cmip-obs","cmip-amip"/)
;compare_type 			= (/"cmip-amip"/)
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
	if(compare_type(icomp) .eq."cmip-amip")then
		fdir 		= (/Dir+"obs/",Dir+"amip/",Dir+"atm/"/)
		if(case(icase).eq."def")then
			casename = (/"","FC5_f09f09_MG15_COSP","BC5_f19g16"/)
		end if
		if(case(icase).eq."new")then
			casename = (/"","new-cldfrac16-1d-MG15-revise-2-COSP","BC5_f19g16_mac2"/)
		end if
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


latS 				= -45
latE 				= 45

lonS				= 180
lonE				= 230



wks   = gsn_open_wks ("x11", "./figure/lat_time_"+var(ivar)+"_"+compare_type(icomp)+"_"+latS+"_"+latE+"_"+lonS+"_"+lonE) 
plot = new((/6/),graphic)


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

;;;;; zonal mean

latS 				= -45
latE 				= 45

lonS				= 180
lonE				= 230


; zonal mean over eastern pacific
var_interp_zm 			= dim_avg_n_Wrap(var_interp(:,:,:,{lonS:lonE}),3)
var_interp_diff_zm 		= dim_avg_n_Wrap(var_interp_diff(:,:,:,{lonS:lonE}),3)

var_interp_zm_rev 		= var_interp_zm(case|:,lat|:,time|:)
var_interp_diff_zm_rev 		= var_interp_diff_zm(case|:,lat|:,time|:)

printVarSummary(var_interp_diff_zm_rev)
printMinMax(var_interp_zm_rev,True)
printMinMax(var_interp_diff_zm_rev,True)


;*************************
; plotting parameters
;*************************
;cmap1 = read_colormap_file("WhiteBlueGreenYellowRed")
;cmap2 = read_colormap_file("ViBlGrWhYeOrRe")
;cmap2 = read_colormap_file("BlueRed")

cmap1				= read_colormap_file("WhiteBlueGreenYellowRed")
cmap2				= read_colormap_file("NCV_blu_red")
;cmap2				= read_colormap_file("NCV_blue_red")
;cmap2				= read_colormap_file("BlueRed")
;cmap2				= read_colormap_file("ViBlGrWhYeOrRe")

plot 				= new(6,graphic)

res                  		= True
res@gsnFrame 			= False
res@gsnDraw 			= False

res@cnFillOn         		= True
res@cnLinesOn        		= False

res@lbOrientation    		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= -0.04

res@vpXF             		= 0.12                ; default is 0.2 (aspect ratio)
res@vpYF             		= 0.8                 ; default is 0.8
res@vpHeightF        		= 0.4                 ; default is 0.6
res@vpWidthF         		= 0.75                ; default is 0.6

res@tmXBMode 			= "Explicit"
res@tmXBValues 			= fspan(0,11,12)
res@tmXBLabels 			= (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

res@lbLabelStride    		= 2

if(lonS.gt.180)then
res@gsnCenterString 		= (360.-lonE)+"W-"+(360-lonS)+"W"
else
res@gsnCenterString 		= lonS+"E-"+(360-lonE)+"W"
end if

res@gsnRightString 			= var(ivar)



res1					= res
res1@cnFillPalette        		= cmap1(:170,:)

res1@cnLevelSelectionMode		= "ExplicitLevels"
res1@cnLevels				= var_levs(ivar,:)


if(compare_type(icomp).eq."cmip-amip")then
	string_cli 			= (/"OBS","AMIP_"+case(icase),"Coupled_"+case(icase)/)
else
	string_cli 			= (/"OBS","DEF","NEW"/)
end if


do icli = 0,2
	res1@gsnLeftString 		= string_cli(icli)
	plot(icli) 			= gsn_csm_lat_time(wks, var_interp_zm_rev(icli,{latS:latE},:), res1) 
end do ; do icli=

;gsn_panel(wks,plot(0:2),(/3,1/),False)

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
	 plot(idiff+3) 			= gsn_csm_lat_time(wks, var_interp_diff_zm_rev(idiff,{latS:latE},:), res2)

end do ; do idff=

; gsn_panel(wks,plot(3:5),(/3,1/),False)


resP 				= True
;resP@gsnPanelLabelBar		= True
resP@amJust			= "TopLeft"
resP@gsnPanelFigureStrings	= (/"a)","b)","c)","d)"/)

gsn_panel(wks,(/plot(0),plot(1),plot(3),plot(5)/),(/2,2/),resP)


resP2				= True
resP2@amJust			= "TopLeft"
resP2@gsnPanelFigureStrings	= (/"a)","b)","c)","d)","e)","f)"/)

gsn_panel(wks,(/plot(0),plot(3),plot(1),plot(4),plot(2),plot(5)/),(/3,2/),resP2)



delete(var_obs)
delete(lat_obs)
delete(lon_obs)

delete(var_def)
delete(lat_def)
delete(lon_def)

delete(var_new)
delete(lat_new)
delete(lon_new)

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
