load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

;%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regrid_diff(var_obs[*][*]:numeric,lat_obs[*]:numeric,lon_obs[*]:numeric,\
                     var_def[*][*]:numeric,lat_def[*]:numeric,lon_def[*]:numeric,\
		     var_new[*][*]:numeric,lat_new[*]:numeric,lon_new[*]:numeric)
begin
;horizontal interpolation
var_interp 			= lat_lon_interpolate_3vars(var_obs,lat_obs,lon_obs, \
						            var_def,lat_def,lon_def, \
				 	                    var_new,lat_new,lon_new)

var_interp_diff 		= var_interp
var_interp_diff(0,:,:) 		= var_interp(1,:,:) - var_interp(0,:,:) ; def-obs
var_interp_diff(1,:,:) 		= var_interp(2,:,:) - var_interp(0,:,:) ; new-obs
var_interp_diff(2,:,:) 		= var_interp(2,:,:) - var_interp(1,:,:) ; new-obs

return([/var_interp,var_interp_diff/])

delete(var_interp)
delete(var_interp_diff)

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin
;*********************************************
; read in data
;*********************************************

season				= (/"ANN"/)
CLDLOW_levs			= (/5,10,15,20,25,30,40,50,60,70,75,80,85,90,95/)
FSNT_levs			= (/25.,50.,75.,100.,125.,150.,175.,200.,225.,250.,275.,300.,325.,350.,375/)
PRECT_levs			= (/.5,1,2,3,4,5,6,7,8,9,10,12,14,15,17/)
PSL_levs			= (/984,988,992,996,1000,1004,1008,1012,1016,1020,1024,1028,1032,1036,1040/)
SST_levs			= (/-1,0,1,3,6,9,12,15,18,20,22,24,26,28,29/)



FSNT_Diff_levs			= (/-80,-60,-40,-30,-20,-10,-5,0,5,10,20,30,40,60,80/)
PSL_Diff_levs			= (/-15,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,15/)
SST_Diff_levs			= (/-5,-4,-3,-2,-1,-.5,-.2,0,.2,.5,1,2,3,4,5/)

PRECT_Diff_levs			 	= (/-5,-4,-3,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,4,5/)
PRECT_Diff2_levs			= (/-4,-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2.,2.5,3.,4/)

CLDLOW_Diff_levs		= (/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)
CLDLOW_Diff2_levs		= (/-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35/)


;var				= (/"PRECT","CLDLOW","SST"/)
;var_levs			= (/PRECT_levs,CLDLOW_levs,SST_levs/)
;varDiff_levs			= (/PRECT_Diff_levs,CLDLOW_Diff_levs,SST_Diff_levs/)

var				= (/"PRECT","CLDLOW"/)
var_levs			= (/PRECT_levs,CLDLOW_levs/)
varDiff_levs			= (/PRECT_Diff_levs,CLDLOW_Diff_levs/)
varDiff2_levs			= (/PRECT_Diff2_levs,CLDLOW_Diff2_levs/)



;var_cmp 			= (/"UV925"/)
var_cmp				= (/"stress"/)

;compare_type 			= (/"amip-obs","cmip-obs","cmip-amip"/)
;compare_type 			= (/"cmip-amip"/)
;compare_type 			= (/"cmip-obs"/)
;compare_type			= (/"amip-obs"/)
compare_type			= (/"cmip-obs"/)

case				= (/"def","new"/)

Dir				= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

do icase=0,dimsizes(case)-1
print(case(icase))

	if(compare_type(icomp) .eq."amip-obs")then
		casename 		= (/"","FC5_f09f09_MG15_COSP","new-cldfrac16-1d-MG15-revise-2-COSP"/)
		fdir 			= (/Dir+"obs/",Dir+"amip/",Dir+"amip/"/)

	end if
	if(compare_type(icomp) .eq."cmip-obs")then
		casename 		= (/"","BC5_f19g16","BC5_f19g16_mac2"/)
		fdir 			= (/Dir+"obs/",Dir+"atm/",Dir+"atm/"/)

	end if
	if(compare_type(icomp) .eq."cmip-amip")then
		fdir 			= (/Dir+"obs/",Dir+"amip/",Dir+"atm/"/)
		if(case(icase).eq."def")then
			casename 	= (/"","FC5_f09f09_MG15_COSP","BC5_f19g16"/)
		end if
		if(case(icase).eq."new")then
			casename 	= (/"","new-cldfrac16-1d-MG15-revise-2-COSP","BC5_f19g16_mac2"/)
		end if
	end if

do is = 0,dimsizes(season)-1

fapend 					= (/"_"+season(is)+"_climo.nc", \
       					    "_"+season(is)+"_climo.nc", \
	   				    "_"+season(is)+"_climo.nc" /)


do ivar=0,dimsizes(var)-1
print(var(ivar))

if(var(ivar).eq."SST")then
	casename(0) 			= "HadISST_PD"
	fdir 				= (/Dir+"obs/", \
	     				    Dir+"ocn/", \
				 	    Dir+"ocn/"/)
end if

if(var(ivar).eq."PRECT")then
	casename(0) 			= "GPCP"
end if

if(var(ivar).eq."CLDLOW".or.var(ivar).eq."CLDMED".or.var(ivar).eq."CLDHGH".or.var(ivar).eq."CLDTOT")then
	casename(0) 			= "CLOUDSAT"
end if

if(var(ivar).eq."FSNS" .or. var(ivar).eq."FLNS" .or. \
   var(ivar).eq."LHFLX" .or. var(ivar).eq."SHFLX" .or. \
   var(ivar).eq."UV925".or.var(ivar).eq."stress".or.\
   var(ivar).eq."PSL".or.\
   var(ivar).eq."NetSfcHeat")then
	casename(0) 			= "JRA25"
end if

if(var(ivar).eq."FSNTOA")then
	casename(0)			= "CERES-EBAF"
end if


;=====================================================================
; For Normal scalar variables
;=====================================================================

print(casename(0))

result_obs				= get_var_2D(fdir(0),fapend(0),casename(0),var(ivar),"obs")
var_obs					= result_obs[0]
lat_obs					= result_obs[1]
lon_obs					= result_obs[2]

result_def				= get_var_2D(fdir(1),fapend(1),casename(1),var(ivar),"model")
var_def					= result_def[0]
lat_def					= result_def[1]
lon_def					= result_def[2]

result_new				= get_var_2D(fdir(2),fapend(2),casename(2),var(ivar),"model")
var_new					= result_new[0]
lat_new					= result_new[1]
lon_new					= result_new[2]

result 					= regrid_diff(var_obs,lat_obs,lon_obs,\
	             				      var_def,lat_def,lon_def,\
		     				      var_new,lat_new,lon_new)
var_interp 				= result[0]
var_interp_diff 			= result[1]

printMinMax(var_interp,True)

delete([/var_obs,lat_obs,lon_obs,var_def,lat_def,lon_def,var_new,lat_new,lon_new/])

;=====================================================================
; For overlaying wind or wind stress field
;=====================================================================

if(var_cmp.eq."UV925".or.var_cmp.eq."stress")then
	if(compare_type(icomp).eq."cmip-obs")then
		casename(0) 			= "JRA25"
		fdir 				= (/Dir+"obs/", \
		         		            Dir+"atm/", \
			 		            Dir+"atm/"/)
	end if	
	if(compare_type(icomp).eq."amip-obs")then
		casename(0) 			= "JRA25"
		fdir 				= (/Dir+"obs/", \
		         		            Dir+"amip/", \
			 		            Dir+"amip/"/)
	end if



end if

if(var_cmp.eq."UV925")then
	result_obs			= get_wind_2D(fdir(0),fapend(0),casename(0),var_cmp,"obs")
	U_obs				= result_obs[0]
	V_obs				= result_obs[1]
	lat_obs				= result_obs[2]
	lon_obs				= result_obs[3]
	
	result_def			= get_wind_2D(fdir(1),fapend(1),casename(1),var_cmp,"model")
	U_def				= result_def[0]
	V_def				= result_def[1]
	lat_def				= result_def[2]
	lon_def				= result_def[3]
	
	result_new			= get_wind_2D(fdir(2),fapend(2),casename(2),var_cmp,"model")
	U_new				= result_new[0]
	V_new				= result_new[1]
	lat_new				= result_new[2]
	lon_new				= result_new[3]
end if


if(var_cmp.eq."stress")then
	result_obs			= get_stress_2D(fdir(0),fapend(0),casename(0),var_cmp,"obs")
	U_obs				= result_obs[0]
	V_obs				= result_obs[1]
	lat_obs				= result_obs[2]
	lon_obs				= result_obs[3]
	
	result_def			= get_stress_2D(fdir(1),fapend(1),casename(1),var_cmp,"model")
	U_def				= result_def[0]
	V_def				= result_def[1]
	lat_def				= result_def[2]
	lon_def				= result_def[3]
	
	result_new			= get_stress_2D(fdir(2),fapend(2),casename(2),var_cmp,"model")
	U_new				= result_new[0]
	V_new				= result_new[1]
	lat_new				= result_new[2]
	lon_new				= result_new[3]
end if

printVarSummary(U_obs)

result 				= regrid_diff(U_obs,lat_obs,lon_obs,\
	             			      U_def,lat_def,lon_def,\
		     			      U_new,lat_new,lon_new)
U_interp 			= result[0]
U_interp_diff 			= result[1]

result 				= regrid_diff(V_obs,lat_obs,lon_obs,\
	             			      V_def,lat_def,lon_def,\
		     			      V_new,lat_new,lon_new)
V_interp 			= result[0]
V_interp_diff 			= result[1]

;U_interp			= U_interp*1.e4 ; change into N/cm2
;V_interp			= V_interp*1.e4 ; change into N/cm2

delete([/U_obs,V_obs,lat_obs,lon_obs,U_def,V_def,lat_def,lon_def,U_new,V_new,lat_new,lon_new/])


;=====================================================================
; Plotting Parameters
;=====================================================================

;latS 			= -70.
;latE 			= 70.
;lonS 			= 140.
;lonE 			= 290.

latS 			= -40.
latE 			= 40.
lonS 			= 140.
lonE 			= 290.



;latS 			= -90.
;latE 			= 90.
;lonS 			= 0.
;lonE 			= 360.

if(var(ivar).eq."PSL".or.var(ivar).eq."SST".or.var(ivar).eq."PRECT")then
	on_vector		= True
else
	on_vector		= False
end if

if(var(ivar).eq."PRECT")then
	on_line			= True
else
	on_line			= False
end if

;=====================================================================
; Setting Resources
;=====================================================================

wks   = gsn_open_wks ("x11", "./figure/lat_lon_"+var(ivar)+"_"+compare_type(icomp)+"_"+case(icase)+"_"+season(is)+"_onvector-"+on_vector+"_online-"+on_line+"_"+latE ) 

plot				= new((/6/),graphic)

;gsn_define_colormap(wks,"BlueDarkRed18")
;cmap1				= read_colormap_file("WhiteBlueGreenYellowRed")
;cmap2				= read_colormap_file("NCV_blu_red")

if(var(ivar).eq."PRECT")then
cmap1				= read_colormap_file("precip_11lev")
cmap2				= read_colormap_file("precip_diff_12lev")
else
cmap1				= read_colormap_file("WhiteBlueGreenYellowRed")
cmap2				= read_colormap_file("NCV_blu_red")
end if




;if(var(ivar).eq."PRECT")then
;
;gsn_define_colormap(wks,"precip_11lev")
;end if


;===============================================
; Plot separate schemes (obs,def,new)
;===============================================


res 				= True
res@gsnFrame 			= False
res@gsnDraw 			= False

res@cnFillOn 			= True
res@cnLinesOn 			= False
res@cnLineLabelsOn 		= False

;----------------------------------------------------------
; Set Map region
;----------------------------------------------------------
res@mpLandFillColor 		= "gray"
res@mpMinLonF 			= lonS
res@mpMaxLonF 			= lonE
res@mpMinLatF 			= latS
res@mpMaxLatF 			= latE
res@mpCenterLonF 		= 180
;----------------------------------------------------------
; Set Labelbar
;----------------------------------------------------------
res@lbOrientation 		= "Vertical"
res@pmLabelBarOrthogonalPosF 	= -0.01
res@lbLabelStride    		= 2                        ; every other label
;----------------------------------------------------------
; Set X and Y axis; tickmarks
;----------------------------------------------------------
res@tmXBTickSpacingF 		= 20.
res@tmXBLabelFontHeightF  	= 0.014
res@tmYLLabelFontHeightF  	= 0.014

;==========================================================
;==========================================================

res1				= res
;res1@cnFillPalette        	= cmap1(:150,:)
res1@cnFillPalette        	= cmap1

;----------------------------------------------------------
; Set range of contour lines
;----------------------------------------------------------
res1@cnLevelSelectionMode        = "ExplicitLevels"
res1@cnLevels                    = var_levs(ivar,:)
;----------------------------------------------------------
; Set center/right/left strings
;----------------------------------------------------------
res1@gsnCenterString 		= season(is)
res1@gsnRightString 		= var(ivar)
if(compare_type(icomp).eq."cmip-amip")then
	string_cli 		= (/"OBS","AMIP_"+case(icase),"Coupled_"+case(icase)/)
else
	string_cli 		= (/"OBS","DEF","NEW"/)
end if

if(on_vector)then
;----------------------------------------------------------
; Set vector resources
;----------------------------------------------------------
res1@gsnScalarContour 		= True

res1@vcRefAnnoOrthogonalPosF 	= -1.0
res1@vcRefLengthF 		= 0.045
res1@vcMinDistanceF		= 0.025
res1@vcLineArrowColor		= "black"           ; change vector color
res1@vcLineArrowThicknessF	= 2.0               ; change vector thickness
res1@vcVectorDrawOrder		= "PreDraw"        ; draw vectors last

;refvec				= max(abs(sqrt(U_interp(:,{latS:latE},{lonS:lonE})^2+V_interp(:,{latS:latE},{lonS:lonE})^2)))/3.
refvec				= 0.08
res1@vcRefMagnitudeF 		= refvec

;===============================================
; Plot separate schemes (obs,def,new)
;===============================================

do icli = 0,2
	res1@gsnLeftString 	= string_cli(icli)
	plot(icli) 		= gsn_csm_vector_scalar_map_ce(wks, U_interp(icli,:,:),V_interp(icli,:,:), var_interp(icli,:,:),res1) 
end do ; do icli=

;gsn_panel(wks,plot(0:2),(/3,1/),False)

else

do icli = 0,2
	res1@gsnLeftString 	= string_cli(icli)
	plot(icli) 		= gsn_csm_contour_map_ce(wks, var_interp(icli,:,:),res1) 
end do ; do icli=

;gsn_panel(wks,plot(0:2),(/3,1/),False)

end if ; on_vector

delete(res1@cnLevels)

;----------------------------------------------------------
; Set center/right/left strings
;----------------------------------------------------------
if(compare_type(icomp).eq."cmip-amip")then
	string_diff = (/"AMIP_"+case(icase)+"-OBS",\
	               "Coupled_"+case(icase)+"-OBS",\
		       "Coupled_"+case(icase)+"-AMIP_"+case(icase)/)
else
	string_diff = (/"DEF-OBS","NEW-OBS","NEW-DEF"/)
end if

;===============================================
; Plot difference plots
;===============================================
res2				= res
res2@cnFillPalette        	= cmap2

;----------------------------------------------------------
; Set vector resources
;----------------------------------------------------------

;refvec				= max(abs(sqrt(U_interp_diff(:,{latS:latE},{lonS:lonE})^2+V_interp_diff(:,{latS:latE},{lonS:lonE})^2)))/3.
;res@vcMinMagnitudeF 		= refvec/3

;if(var(ivar).eq."PRECT")then
;	res@cnFillPalette = "precip_diff_12lev"
;end if

do idiff = 0,2
	res2@gsnLeftString 	= string_diff(idiff)

	; set level values
	res2@cnLevelSelectionMode        = "ExplicitLevels"
	if(idiff.ne.2)then
		res2@cnLevels                    = varDiff_levs(ivar,:)
	else
		res2@cnLevels                    = varDiff2_levs(ivar,:)
	end if

	
	if(on_vector)then
		res2@gsnScalarContour 		= True
		
		res2@vcRefAnnoOrthogonalPosF 	= -1.0
		res2@vcRefLengthF 		= 0.045
		res2@vcMinDistanceF		= 0.025
		res2@vcLineArrowColor		= "black"           ; change vector color
		res2@vcLineArrowThicknessF	= 2.0               ; change vector thickness
		res2@vcVectorDrawOrder		= "PreDraw"        ; draw vectors last

		refvec				= 0.03
		res2@vcRefMagnitudeF 		= refvec

		plot(idiff+3) 		= gsn_csm_vector_scalar_map_ce(wks,U_interp_diff(idiff,:,:),V_interp_diff(idiff,:,:),var_interp_diff(idiff,:,:), res2) 
	else
		plot(idiff+3) 		= gsn_csm_contour_map_ce(wks,var_interp_diff(idiff,:,:), res2) 
	end if	; on_vector

	delete(res2@cnLevels)




end do ; do idff=

;gsn_panel(wks,plot(3:5),(/3,1/),False)
;----------------------------------------------------------
; Panel Plotting
;----------------------------------------------------------
; add one box over eastern pacific that we care more about.
;----------------------------------------------------------

  dum1 = new(4,graphic)
  dum2 = new(4,graphic)

  resp = True
  resp@gsLineThicknessF = 4.0

  resp1	= resp
  resp1@gsLineColor = "white"
  
  resp2	= resp
  resp2@gsLineColor = "grey"
 

  ypts1 = (/-5.0, -5.0, -20.0, -20.0, -5.0/)
  xpts1 = (/180.0, 230.0, 230.0, 180.0, 180.0/)

  ypts2 = (/-5.0, -5.0, -20.0, -20.0, -5.0/)
  xpts2 = (/230.0, 255.0, 255.0, 230.0, 230.0/)

  ; draw each ine separately. Each line must contain two points.

  if(var(ivar).eq."PRECT".or.var(ivar).eq."CLDLOW".or.var(ivar).eq."SST")then
  do i = 0,3
  	dum1(i) = gsn_add_polyline(wks,plot(5),xpts1(i:i+1),ypts1(i:i+1),resp1)
  	dum2(i) = gsn_add_polyline(wks,plot(5),xpts2(i:i+1),ypts2(i:i+1),resp2)
  end do
  end if

;----------------------------------------------------------


;gsn_panel(wks,plot,(/2,3/),False)

resP 				= True
;resP@gsnPanelLabelBar		= True
resP@amJust			= "TopLeft"
resP@gsnPanelFigureStrings	= (/"a)","b)","c)","d)"/)

;gsn_panel(wks,(/plot(0),plot(1),plot(3),plot(5)/),(/2,2/),resP)

resP2				= True
resP2@amJust			= "TopLeft"
resP2@gsnPanelFigureStrings	= (/"a)","b)","c)","d)","e)","f)"/)

gsn_panel(wks,(/plot(0),plot(3),plot(1),plot(4),plot(2),plot(5)/),(/3,2/),resP2)

;delete(res)

delete(res1)
delete(res2)


delete(var_interp)
delete(var_interp_diff)

delete(plot)

delete(cmap1)
delete(cmap2)

end do; do ivar=

end do ; do is=

if(compare_type(icomp) .eq."amip-obs" .or. compare_type(icomp) .eq."cmip-obs")then
	break
end if

end do; do icase=

end do; do icomp=
                               
end
