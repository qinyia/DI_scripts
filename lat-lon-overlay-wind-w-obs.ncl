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

;var				=(/"PRECT","CLDTOT","CLDHGH","CLDMED","CLDLOW"/)
;var				=(/"FSNS","FLNS","LHFLX","SHFLX"/)
;var				=(/"FSNS"/)
;var				=(/"LHFLX"/)
;var				= (/"PRECT"/)
;var				= (/"NetSfcHeat"/)
var				= (/"CLDLOW","FSNT"/)

;var				=(/"PRECT"/)
;var				= (/"SST","CLDLOW","FSNS","NetSfcHeat","PSL"/)
;var				= (/"PSL"/)
;var				= (/"SST"/)

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


;=====================================================================
; For Normal scalar variables
;=====================================================================

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

latS 			= -40.
latE 			= 40.
lonS 			= 140.
lonE 			= 290.

;latS 			= -90.
;latE 			= 90.
;lonS 			= 0.
;lonE 			= 360.

if(var(ivar).eq."PSL".or.var(ivar).eq."SST")then
	on_vector		= True
else
	on_vector		= False
end if

;=====================================================================
; Setting Resources
;=====================================================================

wks   = gsn_open_wks ("x11", "./figure/lat_lon_"+var(ivar)+"_"+compare_type(icomp)+"_"+case(icase)+"_"+season(is)+"_onvector-"+on_vector ) 

plot				= new((/6/),graphic)



;===============================================
; Plot separate schemes (obs,def,new)
;===============================================


res 				= True
res@gsnFrame 			= False
res@gsnDraw 			= False
res@cnFillPalette 		= "BlueDarkRed18"

if(var(ivar).eq."PRECT")then
res@cnFillPalette 		= "precip4_11lev"
end if

;res@gsnSpreadColors 		= True

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
;----------------------------------------------------------
; Set Grid lines
;----------------------------------------------------------
;res@mpGridAndLimbOn        	= True               ; turn on grid lines
;res@mpGridLineDashPattern 	= 5
;res@mpGridLatSpacingF 		= 10.
;res@mpGridLonSpacingF 		= 20.
;res@mpGridMaskMode 		= "MaskNotOcean"
;----------------------------------------------------------
; Set range of contour lines
;----------------------------------------------------------
maxval                          = max(var_interp(:,{latS:latE},{lonS:lonE}))
minval                          = min(var_interp(:,{latS:latE},{lonS:lonE}))

spcval                          = (maxval-minval)/11

data_levels                     = fspan(minval,maxval,11)

res@cnLevelSelectionMode        = "ExplicitLevels"
res@cnLevels                    = data_levels

res@cnLineLabelFormat           = "0@;*.3f"                      ; Only one value after decimal point
if(abs(spcval).lt.1.e-4)then
        res@lbLabelStrings	= sprintf("%5.1g",data_levels)   ; Format the labelbar labels
else
        res@lbLabelStrings	= sprintf("%5.1f",data_levels)   ; Format the labelbar labels
end if
;----------------------------------------------------------
; Set center/right/left strings
;----------------------------------------------------------
res@gsnCenterString = season(is)
res@gsnRightString = var(ivar)
if(compare_type(icomp).eq."cmip-amip")then
	string_cli 		= (/"OBS","AMIP_"+case(icase),"Coupled_"+case(icase)/)
else
	string_cli 		= (/"OBS","DEF","NEW"/)
end if

if(on_vector)then
;----------------------------------------------------------
; Set vector resources
;----------------------------------------------------------
res@gsnScalarContour 		= True

res@vcRefAnnoOrthogonalPosF 	= -1.0
res@vcRefLengthF 		= 0.045
res@vcMinDistanceF		= 0.025
res@vcLineArrowColor		= "black"           ; change vector color
res@vcLineArrowThicknessF	= 2.0               ; change vector thickness
res@vcVectorDrawOrder		= "PostDraw"        ; draw vectors last

;refvec				= max(abs(sqrt(U_interp(:,{latS:latE},{lonS:lonE})^2+V_interp(:,{latS:latE},{lonS:lonE})^2)))/3.
refvec				= 0.08
res@vcRefMagnitudeF 		= refvec

;===============================================
; Plot separate schemes (obs,def,new)
;===============================================

do icli = 0,2
	res@gsnLeftString 	= string_cli(icli)
	plot(icli) 		= gsn_csm_vector_scalar_map_ce(wks, U_interp(icli,:,:),V_interp(icli,:,:), var_interp(icli,:,:),res) 
end do ; do icli=

gsn_panel(wks,plot(0:2),(/3,1/),False)

else

do icli = 0,2
	res@gsnLeftString 	= string_cli(icli)
	plot(icli) 		= gsn_csm_contour_map_ce(wks, var_interp(icli,:,:),res) 
end do ; do icli=

;gsn_panel(wks,plot(0:2),(/3,1/),False)

end if ; on_vector

delete(maxval)
delete(minval)
delete(spcval)
delete(data_levels)
delete(res@cnLevels)

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
;----------------------------------------------------------
; Set vector resources
;----------------------------------------------------------

;refvec				= max(abs(sqrt(U_interp_diff(:,{latS:latE},{lonS:lonE})^2+V_interp_diff(:,{latS:latE},{lonS:lonE})^2)))/3.
;res@vcMinMagnitudeF 		= refvec/3

if(var(ivar).eq."PRECT")then
	res@cnFillPalette = "precip_diff_12lev"
end if

do idiff = 0,2
	res@gsnLeftString 	= string_diff(idiff)

	if(idiff.ne.2)then
		;----------------------------------------------------------
		; Set range of contour lines
		;----------------------------------------------------------
		maxval                          = max(var_interp_diff(0:1,{latS:latE},{lonS:lonE}))
		minval                          = min(var_interp_diff(0:1,{latS:latE},{lonS:lonE}))
	else
		maxval                          = max(var_interp_diff(2,{latS:latE},{lonS:lonE}))
		minval                          = min(var_interp_diff(2,{latS:latE},{lonS:lonE}))
	end if

	spcval                          = (maxval-minval)/11.
	
	tmp                             = (/abs(maxval),abs(minval)/)
	exmval                          = 1./2*(max(tmp)+min(tmp))
	data_levels                     = fspan(-1*exmval,exmval,11)
	
	res@cnLevelSelectionMode        = "ExplicitLevels"
	res@cnLevels                    = data_levels
	
	res@cnLineLabelFormat           = "0@;*.3f"                      ; Only one value after decimal point
	if(abs(spcval).lt.1.e-4)then
	        res@lbLabelStrings	= sprintf("%5.1g",data_levels)   ; Format the labelbar labels
	else
	        res@lbLabelStrings	= sprintf("%5.1f",data_levels)   ; Format the labelbar labels
	end if

	if(on_vector)then
		refvec				= 0.03
		res@vcRefMagnitudeF 		= refvec

		plot(idiff+3) 		= gsn_csm_vector_scalar_map_ce(wks,U_interp_diff(idiff,:,:),V_interp_diff(idiff,:,:),var_interp_diff(idiff,:,:), res) 
	else
		plot(idiff+3) 		= gsn_csm_contour_map_ce(wks,var_interp_diff(idiff,:,:), res) 
	end if	; on_vector

	delete(maxval)
	delete(minval)
	delete(spcval)
	delete(data_levels)
	delete(res@cnLevels)




end do ; do idff=

;gsn_panel(wks,plot(3:5),(/3,1/),False)
;----------------------------------------------------------
; Panel Plotting
;----------------------------------------------------------
; add one box over eastern pacific that we care more about.
;----------------------------------------------------------

  dum1 = new(4,graphic)

  resp = True
  resp@gsLineColor = "white"
  resp@gsLineThicknessF = 4.0

  ypts = (/-5.0, -5.0, -20.0, -20.0, -5.0/)
;  xpts = (/180.0, 220.0, 220.0, 180.0, 180.0/)
  xpts = (/180.0, 230.0, 230.0, 180.0, 180.0/)

  ; draw each ine separately. Each line must contain two points.

  if(var(ivar).eq."PRECT")then
  do i = 0,3
  	dum1(i) = gsn_add_polyline(wks,plot(5),xpts(i:i+1),ypts(i:i+1),resp)
  end do
  end if

;----------------------------------------------------------


;gsn_panel(wks,plot,(/2,3/),False)

resP 				= True
;resP@gsnPanelLabelBar		= True
resP@amJust			= "TopLeft"
resP@gsnPanelFigureStrings	= (/"a)","b)","c)","d)"/)

gsn_panel(wks,(/plot(0),plot(1),plot(3),plot(5)/),(/2,2/),resP)

;gsn_panel(wks,(/plot(0),plot(3),plot(1),plot(4),plot(2),plot(5)/),(/3,2/),resP)

delete(res)


delete(var_interp)
delete(var_interp_diff)

delete(plot)

end do; do ivar=

end do ; do is=

if(compare_type(icomp) .eq."amip-obs" .or. compare_type(icomp) .eq."cmip-obs")then
	break
end if

end do; do icase=

end do; do icomp=
                               
end
