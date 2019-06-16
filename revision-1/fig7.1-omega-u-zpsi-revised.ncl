load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

load "./functions_lat_interpolate.ncl"

function zpsi(lev[*]:numeric,lat[*]:numeric,lon[*]:numeric,\
              U_1[*][*][*]:numeric,V_1[*][*][*]:numeric,OMEGA_1[*][*][*],\
	      gw[*]:numeric,latS:float,latE:float)

begin
;==================================================================
; set parameters...
;==================================================================
g			= 9.81		; gravity
a			= 6.37*10^6	; radius of earth
pi		 	= 3.1428	; pi

nlev			= dimsizes(lev)
nlat			= dimsizes(lat)
nlon			= dimsizes(lon)

;==================================================================
; get average over specific zonal ranges with raw wind: u and omega
;==================================================================

OMEGA_1_mm 		= dim_avg_wgt_n_Wrap(OMEGA_1(:,{latS:latE},:),gw({latS:latE}),0,1)
OMEGA_1_mm&lev@units	= "hPa"
printVarSummary(OMEGA_1_mm)
printMinMax(OMEGA_1_mm,True)

U_1_mm 			= dim_avg_wgt_n_Wrap(U_1(:,{latS:latE},:),gw({latS:latE}),0,1)
U_1_mm&lev@units	= "hPa"
printVarSummary(U_1_mm)
printMinMax(U_1_mm,True)

;==================================================================
; get velocity potential: u=d()/dx; v = d()/dy
; and divergent wind: u & v
;==================================================================
sfvp_1 			= new((/2,nlev,nlat,nlon/),float)

sfvp_1(0,:,:,:)		= U_1
sfvp_1 			= uv2sfvpG(U_1,V_1)
printVarSummary(sfvp_1)
;printMinMax(sfvp_1,True)

; sfvp(0,...) is stream function
; sfvp(1,...) is velocity potential
vp_1 			= sfvp_1(1,:,:,:)

DvpDx_1 		= new((/nlev,nlat,nlon/),float)
DvpDx_1			= vp_1
DvpDx_1 		= doubletofloat(center_finite_diff_n(vp_1,lon*111000,False,0,2))

;printVarSummary(DvpDx_1)
;printMinMax(DvpDx_1,True)

DvpDx_1_mm 		= dim_avg_wgt_n_Wrap(DvpDx_1(:,{latS:latE},:),gw({latS:latE}),0,1)
DvpDx_1_mm&lev@units	= "hPa"
;printVarSummary(DvpDx_1_mm)

;==================================================================
; calculate zonal stream function: Walker Circulation
;==================================================================
dpsi 			= new(  (/dimsizes(lev),dimsizes(lon)/) , "float")
psi 			= new(  (/dimsizes(lev),dimsizes(lon)/) , "double")
dpsidphi		= new(  (/dimsizes(lev),dimsizes(lon)/) , "double")
omega			= new(  (/dimsizes(lev),dimsizes(lon)/) , "double")

; initialize these arrays at zeros
psi 			= 0.
omega 			= 0.

;since all of these values use the same level/latitude coordinates
;you can copy v to each other variable
copy_VarCoords(DvpDx_1_mm,dpsi)
copy_VarCoords(DvpDx_1_mm,psi)
copy_VarCoords(DvpDx_1_mm,dpsidphi)
copy_VarCoords(DvpDx_1_mm,omega)

;printVarSummary(dpsi)

; get dpsi
do ilev=0,dimsizes(lev)-1
	do ilon = 0,dimsizes(lon)-1
		dpsi(ilev,ilon)	= DvpDx_1_mm(ilev,ilon)*2*pi*a/g
	end do
end do

do ilev=1,dimsizes(lev)-1
	do ilon=0,dimsizes(lon)-1
		psi(ilev,ilon) = psi(ilev-1,ilon)+dpsi(ilev,ilon)*((lev(ilev)-lev(ilev-1))*100)

	end do
end do

;printVarSummary(dpsi)
;printMinMax(dpsi,True)
;printVarSummary(psi)
;printMinMax(psi,True)

; psi is calculated!
; the next step is to calculate omega which we will do a forward and backwards differencing
; for the first and last values, but centerred differencing for
; all the other calculations

; we do this forwards/bacwards/ and centered differencing on the d(psi)/d(phi) term

do ilev=0,dimsizes(lev)-1
	dpsidphi(ilev,0)		= (psi(ilev,1)-psi(ilev,0))/(lon(1)-lon(0))*(180/pi)
	dpsidphi(ilev,nlon-1)		= (psi(ilev,nlon-1)-psi(ilev,nlon-2))/(lon(1)-lon(0))*(180/pi)

	do ilon=1,nlon-2
		dpsidphi(ilev,ilon)	= (psi(ilev,ilon+1)-psi(ilev,ilon-1))/(lon(1)-lon(0))*(180/pi)
	end do
end do

;printVarSummary(dpsidphi)
;printMinMax(dpsidphi,True)

; Finally run one last set of do lops to calculate omega
do ilev=1,dimsizes(lev)-1
	do ilon=0,dimsizes(lon)-1
		omega(ilev,ilon)	= -g/(2*pi*a^2)*dpsidphi(ilev,ilon)
	end do
end do

;printVarSummary(omega)
;printMinMax(omega,True)

; This step is just to adjust the values so they can easily be plotted. Proper units are cited on the figure.
psi 					= psi*10^-10
;omega					= omega*10^3
omega					= omega*864

;printMinMax(psi,True)
;printMinMax(omega,True)

return([/psi,omega,OMEGA_1_mm,U_1_mm/])

end

;==============================================================
; Main programing starts here....
;==============================================================

begin

;***************file 1********************
;fdir 			= "/home/lyl/WORK3/qinyi/double-ITCZ/data/amip/"
;case 			= (/"FC5_f19f19_MG15_cosp","FC5_f19f19_MG15_mac3_cosp"/)

season			= (/"ANN"/)

;compare_type 		= (/"amip-obs","cmip-obs"/)
compare_type 		= (/"cmip-obs"/)


endyear			= (/"0004","0006","0008","0012"/)

Dir			= "/home/lyl/WORK3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

sstvt_case		= (/"BC5_f19g16_mac2_cosp","BC5_f19g16_sstvt_0.52",\
			    "BC5_f19g16_sstvt_enlarge_region_0.75","BC5_f19g16_sstvt_abs_100percent","BC5_f19g16_sstvt_0.85"/)
;sstvt_case		= (/"BC5_f19g16_sstvt_0.85"/)
;sstvt_case		= (/"BC5_f19g16_sstvt_0.52"/)

do isstcase = 0,dimsizes(sstvt_case)-1
if(compare_type(icomp) .eq."cmip-obs")then
	fdir 		= (/Dir+"obs/",Dir+"atm/",Dir+"atm/"/)
	case		= (/"ERAI","BC5_f19g16_cosp",sstvt_case(isstcase)/)
end if

if(compare_type(icomp) .eq."amip-obs")then
	case		= (/"ERAI","FC5_f19f19_MG15_cosp","FC5_f19f19_MG15_mac3_cosp"/)
	fdir 		= (/Dir+"obs/",Dir+"amip/",Dir+"amip/"/)
end if

do is			= 0, dimsizes(season)-1

do iend 		= 0, dimsizes(endyear)-1 

f0 			= addfile(fdir(0)+case(0)+"_"+season(is)+"_climo.nc","r")
f1 			= addfile(fdir(1)+case(1)+"_"+season(is)+"_climo_"+endyear(iend)+".nc","r")
f2 			= addfile(fdir(2)+case(2)+"_"+season(is)+"_climo_"+endyear(iend)+".nc","r")

lev_obs			= f0->lev
lat_obs			= f0->lat
lon_obs			= f0->lon
U_obs 			= f0->U(0,:,:,:)
V_obs 			= f0->V(0,:,:,:)
OMEGA_obs		= f0->OMEGA(0,:,:,:) ; mb/day
OMEGA_obs		= OMEGA_obs/864. ; Pa/s

lev_def			= f1->lev
lat_def			= f1->lat
lon_def			= f1->lon
U_def 			= f1->U(0,:,:,:)
V_def 			= f1->V(0,:,:,:)
OMEGA_def 		= f1->OMEGA(0,:,:,:)

lev_new			= f2->lev
lat_new			= f2->lat
lon_new			= f2->lon
U_new 			= f2->U(0,:,:,:)
V_new 			= f2->V(0,:,:,:)
OMEGA_new 		= f2->OMEGA(0,:,:,:)


; regrid all data to the same grid information
lat_spec		= lat_def;fspan(-90,90,72)
lon_spec		= lon_def;fspan(0,360,144)
lev_spec		= lev_def
; be careful: the streamfunction is sensitive to the number of vertical levels.
; here I use the model lev as the aimed level rather than the specified one.

U_regrid		= lev_lat_lon_interpolate_3vars_specgrid(U_obs,lat_obs,lon_obs,lev_obs,\
							         U_def,lat_def,lon_def,lev_def,\
							         U_new,lat_def,lon_def,lev_new,\
							         lat_spec,lon_spec,lev_spec)

V_regrid		= lev_lat_lon_interpolate_3vars_specgrid(V_obs,lat_obs,lon_obs,lev_obs,\
							         V_def,lat_def,lon_def,lev_def,\
							         V_new,lat_def,lon_def,lev_new,\
							         lat_spec,lon_spec,lev_spec)

OMEGA_regrid		= lev_lat_lon_interpolate_3vars_specgrid(OMEGA_obs,lat_obs,lon_obs,lev_obs,\
							         OMEGA_def,lat_def,lon_def,lev_def,\
							         OMEGA_new,lat_def,lon_def,lev_new,\
							         lat_spec,lon_spec,lev_spec)

printVarSummary(OMEGA_regrid)

U_1			= U_regrid(0,:,:,:)
V_1			= V_regrid(0,:,:,:)
OMEGA_1			= OMEGA_regrid(0,:,:,:)

U_2			= U_regrid(1,:,:,:)
V_2			= V_regrid(1,:,:,:)
OMEGA_2			= OMEGA_regrid(1,:,:,:)

U_3			= U_regrid(2,:,:,:)
V_3			= V_regrid(2,:,:,:)
OMEGA_3			= OMEGA_regrid(2,:,:,:)

;1111111111111111111111111111111

latS 			= -20.
latE 			= -5.
lonS			= 180.
lonE			= 260.

g			= 9.81		; gravity
a			= 6.37*10^6	; radius of earth
pi		 	= 3.1428	; pi

lat_here		= U_1&lat
rad			= 4.0*atan(1.0)/180.
clat			= cos(lat_here*rad)
clat!0			= "lat"
clat&lat		= lat_here

lev			= U_1&lev
lat			= U_1&lat
lon			= U_1&lon

;==================================================================
; get zonal stream function from subroutine
;==================================================================
result 			= zpsi(lev, lat, lon, U_1, V_1, OMEGA_1, clat, latS, latE)
psi_1			= result[0]
omega_1			= result[1]
OMEGA_1_mm		= result[2]
U_1_mm			= result[3]

result 			= zpsi(lev, lat, lon, U_2, V_2, OMEGA_2, clat, latS, latE)
psi_2			= result[0]
omega_2			= result[1]
OMEGA_2_mm		= result[2]
U_2_mm			= result[3]

result 			= zpsi(lev, lat, lon, U_3, V_3, OMEGA_3, clat, latS, latE)
psi_3			= result[0]
omega_3			= result[1]
OMEGA_3_mm		= result[2]
U_3_mm			= result[3]


printMinMax(psi_1,True)
printMinMax(psi_2,True)
printMinMax(psi_3,True)

;================================================================
; Plot starts here...
;================================================================
res 					= True

res@gsnDraw 				= False
res@gsnFrame 				= False

res@trYMinF 				= 100.
res@tiYAxisString			= "Pressure (hPa)"
res@tiYAxisFontHeightF			= 0.0325

res@trXMinF				= lonS
res@trXMaxF				= lonE

res@vpWidthF				= 0.7
res@vpHeightF				= 0.4

res@gsnPresHgtHeightLabelOn		= False
res@tmYROn				= False

pmslres 				= res 
lrhes 					= res
psires 					= res
omegres 				= res

;---------------------------------------------------
;---------------------------------------------------
; This first one is for streamfunction - (psi)
;---------------------------------------------------
;---------------------------------------------------
psires 					= True

psires@gsnLeftString 			= ""
psires@gsnRightString 			= ""
psires@gsnStringFontHeightF		= 0.0325

psires@pmTickMarkDisplayMode 		= "Always"
;---------------------------------------------------
; shaded or lines
;---------------------------------------------------
psires@cnFillOn 			= False
psires@cnLinesOn 			= True                       
psires@cnInfoLabelOn         		= False          
psires@cnLineLabelsOn 			= True
;psires@cnFillDrawOrder      		= "PostDraw"

;---------------------------------------------------
; title
;---------------------------------------------------
psires@tiMainFontHeightF 		= 0.012              
;---------------------------------------------------
; label bar 
;---------------------------------------------------
psires@lbLabelAutoStride 		= True
psires@lbLabelOffsetF 			= 0.1
;psires@lbBoxLinesOn          		= False
;psires@lbTitleString 			= "[10~S~-5~N~s~S~-1~N~]"
;psires@lbTitlePosition 		= "Right"
;psires@lbTitleDirection 		= "Across"
;---------------------------------------------------
; line feature: thickness, color, et al
;---------------------------------------------------
psires@cnLineThicknessF 		= 1.5
psires@gsnContourNegLineDashPattern 	= 2
;psires@gsnContourZeroLineThicknessF 	= 6.0
; hiding the zero line
psires@gsnContourZeroLineThicknessF	= 0.0

psires@tmXBOn				= False
psires@tmYLOn				= False
psires@tmYROn				= False

;---------------------------------------------------
; line range
;---------------------------------------------------
psires@cnLevelSelectionMode 		= "ManualLevels"
psires@cnMinLevelValF       		= -15.
psires@cnMaxLevelValF       		= 15. 
psires@cnLevelSpacingF      		= 3. 

psires@cnLineLabelPlacementMode		= "Computed"
psires@cnLineLabelFontHeightF		= 0.03

;psires@cnExplicitLineLabelsOn		= True
;psires@cnLineLabelStrings		= (/6,9,12/)

;---------------------------------------------------
;---------------------------------------------------
; This first one is for omega 
;---------------------------------------------------
;---------------------------------------------------

omegres 				= True

omegres@gsnLeftString 			= ""
omegres@gsnRightString 			= ""

omegres@pmTickMarkDisplayMode 		= "Always"

omegres@cnFillOn 			= True
omegres@cnLinesOn 			= False
omegres@cnInfoLabelOn         		= False         
omegres@cnLineLabelsOn 			= False                
omegres@tiMainFontHeightF 		= 0.012             
omegres@cnFillDrawOrder      		= "PreDraw"
;---------------------------------------------------
; label bar 
;---------------------------------------------------
omegres@lbLabelAutoStride 		= True
;omegres@lbBoxLinesOn          		= False
;omegres@lbTitleString 			= "Omega [10~S~-3~N~ hPa s~S~-1~N~]"
;omegres@lbTitlePosition 		= "Right"
;omegres@lbTitleDirection 		= "Across"

omegres@lbTitleString			= "hPa/day"
;omegres@lbTitleJust			= "BottomRight"

;omegres@lbTitlePosition		= "Right"
;omegres@lbTitleDirection		= "Across"

omegres@lbTitleFontHeightF		= 0.03
;---------------------------------------------------
; line range
;---------------------------------------------------
omegres@cnLevelSelectionMode 		= "ManualLevels"
omegres@cnMinLevelValF       		= -50.
omegres@cnMaxLevelValF       		= 50.
omegres@cnLevelSpacingF      		= 5. 

;*******************((CREATE PLOTS))******************
;*****************************************************
wks_type 				= "pdf"
;wks_type@wkWidth 			= 1250
;wks_type@wkHeight 			= 1250

plotName 				= "fig7.1-omega-u-zpsi_"+compare_type(icomp)+"_"+season(is)+"_"+latS+"to"+latE+"_revised_w_obs_"+case(2)+"_"+endyear(iend)
wks_name 				= "./figure/"+plotName
wks 					= gsn_open_wks(wks_type,wks_name)


gsn_define_colormap(wks,"BlueYellowRed")
;gsn_define_colormap(wks,"BlueRed")

omegres@gsnRightStringFontHeightF 	= .025
omegres@gsnRightStringOrthogonalPosF 	= .02

omegres@tmXBLabelFontHeightF		= 0.0325
omegres@tmYLLabelFontHeightF		= 0.0325
omegres@tmYRLabelFontHeightF		= 0.0325
omegres@gsnStringFontHeightF		= 0.0325
omegres@lbLabelFontHeightF		= 0.0325
omegres@tmYLMode                    	= "Explicit"
omegres@tmYLValues                  	= (/850,700,500,300,200/)
omegres@tmYLLabels                  	= (/850,700,500,300,200/)

omegres@tmXBTickSpacingF		= 15

plot					= new(3,graphic)

plot(0)  				= gsn_csm_pres_hgt(wks,psi_1,psires)
plot(1)  				= gsn_csm_pres_hgt(wks,psi_2,psires)
plot(2)					= gsn_csm_pres_hgt(wks,psi_3,psires)

omegres@gsnLeftString			= "OBS"
plot_ov1  				= gsn_csm_pres_hgt(wks,omega_1,omegres)
omegres@gsnLeftString			= "DEF"
plot_ov2  				= gsn_csm_pres_hgt(wks,omega_2,omegres)
omegres@gsnLeftString			= "NEW"
plot_ov3  				= gsn_csm_pres_hgt(wks,omega_3,omegres)

overlay(plot(0),plot_ov1)
overlay(plot(1),plot_ov2)
overlay(plot(2),plot_ov3)


;gsn_panel(wks,plot,(/1,3/),False)

;---------------------------------------------------
; line range
;---------------------------------------------------
if(compare_type(icomp).eq."amip-obs")then
	psires@cnLevelSelectionMode 		= "ManualLevels"
	psires@cnMinLevelValF       		= -4.
	psires@cnMaxLevelValF       		= 4.
	psires@cnLevelSpacingF      		= 1.

	omegres@cnLevelSelectionMode 		= "ManualLevels"
	omegres@cnMinLevelValF       		= -15.
	omegres@cnMaxLevelValF       		= 15.
	omegres@cnLevelSpacingF      		= 5.
end if

if(compare_type(icomp).eq."cmip-obs")then
	psires@cnLevelSelectionMode 		= "ManualLevels"
	psires@cnMinLevelValF       		= -10.
	psires@cnMaxLevelValF       		= 10.
	psires@cnLevelSpacingF      		= 1.

	omegres@cnLevelSelectionMode 		= "ManualLevels"
	omegres@cnMinLevelValF       		= -20.
	omegres@cnMaxLevelValF       		= 20.
	omegres@cnLevelSpacingF      		= 5.
end if


omega_21				= omega_1
omega_21				= (/omega_2-omega_1/)
psi_21					= psi_1
psi_21					= (/psi_2-psi_1/)

omega_31				= omega_1
omega_31				= (/omega_3-omega_1/)
psi_31					= psi_1
psi_31					= (/psi_3-psi_1/)

omega_32				= omega_2
omega_32				= (/omega_3-omega_2/)
psi_32					= psi_2
psi_32					= (/psi_3-psi_2/)

plotd					= new(3,graphic)


omegres@gsnLeftString			= "DEF-OBS"
plotd(0) 				= gsn_csm_pres_hgt(wks,psi_21,psires)
plot_ov1  				= gsn_csm_pres_hgt(wks,omega_21,omegres)
overlay(plotd(0),plot_ov1)

omegres@gsnLeftString			= "NEW-OBS"
plotd(1) 				= gsn_csm_pres_hgt(wks,psi_31,psires)
plot_ov2  				= gsn_csm_pres_hgt(wks,omega_31,omegres)
overlay(plotd(1),plot_ov2)

omegres@gsnLeftString			= "NEW-DEF"
plotd(2) 				= gsn_csm_pres_hgt(wks,psi_32,psires)
plot_ov3  				= gsn_csm_pres_hgt(wks,omega_32,omegres)
overlay(plotd(2),plot_ov3)


resP				= True
resP@amJust                     = "BottomLeft"
;resP@gsnPanelFigureStrings      = (/"a)","b)","c)","d)","e)","f)"/)
resP@gsnPanelXWhiteSpacePercent	= 5.
resP@gsnPanelYWhiteSpacePercent	= 5.


gsn_panel(wks,ndtooned((/plot,plotd/)),(/2,3/),resP)


print( "Figures are finished")
system( "date" )

end do ; do iend=

end do ; do is=

end do ; do isstcase=

end do ; do icomp=
end

