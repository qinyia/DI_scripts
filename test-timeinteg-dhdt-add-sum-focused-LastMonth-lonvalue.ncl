load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************
season				=(/"ANN"/)

;DoTimeInteg			= True
DoTimeInteg			= False

varTimeDepthLatLon		=(/"TEMP","UVEL","VVEL","WVEL"/)
varTimeLatLon			=(/"HMXL","SHF"/)

varTerms_tmp			= (/"netFlux","u_dTdx","v_dTdy","hordiff","w_dTdz","vertdiff3","Qpenflux"/)
varTermsDiff_tmp		= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","hordiff_diff","w_dTdz_diff","vertdiff3_diff","Qpenflux_diff"/)

varTerms			= (/"netFlux","u_dTdx","v_dTdy","hordiff","w_dTdz","vertdiff3"/)
varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","hordiff_diff","w_dTdz_diff","vertdiff3_diff"/)
colors				= (/"brown4","aquamarine","cadetblue","blue","magenta2","purple2"/)

dashes				= fspan(0,0,dimsizes(varTerms))
thick				= fspan(3,3,dimsizes(varTerms))

; add new terms
addterms		= 7
legends_all		=  array_append_record(varTerms,(/"h_dhdt",\
							"sumof4","no h_dhdt","no_Qpenflux","no_vertdiff",\
							"dSSTdt","residual"/),0)

colors_all		=  array_append_record(colors,(/"green",\
							"grey","burlywood3","bisque2","antiquewhite1",\
							"black","grey"/),0)

dashes_all		=  array_append_record(dashes,(/0.,\
							0.,0.,0.,0.,\
							0.,14./),0)
thick_all		=  array_append_record(thick,(/3.,\
						       4.,3.,3.,3.,\
						       4.,4./),0)
;==================================================================
; select variables that we need to draw plots
;(/"netFlux","u_dTdx","v_dTdy","hordiff","w_dTdz","vertdiff3","h_dhdt",\
;  "sumof4","no h_dhdt","no Qpenflux", "no vertdiff","dSSTdt","residual"/)
; (/0,1,2,3,4,5,6,\
;   7,8,9,10,11,12/)
;==================================================================
;vector				=	(/0,1,2,3,4,5,6,7,11,12/)
;vector1				=	(/0,1,2,3,4/)
;vector2				=	(/5,6,7,11,12/)
;vector_empha			=	(/7,11,12/)

vector				=	(/0,1,2,7,11/)
vector1				=	(/0,1,2/)
vector2				=	(/7,11/)
vector_empha			=	(/7,11/)

;==================================================================
;==================================================================



do is = 0,dimsizes(season)-1

;-----------------------------------------------
; Do time integrals...
;-----------------------------------------------
if(DoTimeInteg)then

;========================================================
; read NC file
;========================================================
ncdf = addfile("test-cal-terms-vertdiff3-Qpenflux.nc" ,"r")

dSizes                          = getfiledimsizes(ncdf)  

; directly read all variables
vars_Terms_tmp         		= new((/dimsizes(varTerms_tmp),dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_Terms_tmp!0		= "vars"
vars_Terms_tmp&vars		= ispan(0,dimsizes(varTerms_tmp)-1,1)

vars_TermsDiff_tmp         	= new((/dimsizes(varTermsDiff_tmp),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_TermsDiff_tmp!0		= "vars"
vars_TermsDiff_tmp&vars		= ispan(0,dimsizes(varTermsDiff_tmp)-1,1)

do ivar=0,dimsizes(varTerms_tmp)-1
        vars_Terms_tmp(ivar,:,:,:,:)           =ncdf->$varTerms_tmp(ivar)$ 
end do

do ivar=0,dimsizes(varTermsDiff_tmp)-1
	vars_TermsDiff_tmp(ivar,:,:,:)		= ncdf->$varTermsDiff_tmp(ivar)$
end do

; save all read variables and other calculated variables
vars_Terms         	= new((/dimsizes(varTerms)+addterms,dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_Terms!0		= "vars"
vars_Terms&vars		= ispan(0,dimsizes(varTerms)-1+addterms,1)

vars_TermsDiff         	= new((/dimsizes(varTermsDiff)+addterms,dSizes(1),dSizes(2),dSizes(3)/),float)
vars_TermsDiff!0	= "vars"
vars_TermsDiff&vars	= ispan(0,dimsizes(varTermsDiff)-1+addterms,1)

; put all necessary variables in new arrays: ignore Qpen_flux and replace netFlux with updated net Flux
vars_Terms(:dimsizes(varTerms)-1,:,:,:,:)		= vars_Terms_tmp(:dimsizes(varTerms_tmp)-2,:,:,:,:)
vars_TermsDiff(:dimsizes(varTerms)-1,:,:,:)		= vars_TermsDiff_tmp(:dimsizes(varTerms_tmp)-2,:,:,:)
printVarSummary(vars_Terms)

dims				= dimsizes(vars_Terms)
nvars				= dims(0)
ncase				= dims(1)
ntime				= dims(2)
nlat				= dims(3)
nlon				= dims(4)

; Get the total input heat flux
Q_flux			= vars_Terms_tmp(0,:,:,:,:)
Qpen_flux		= vars_Terms_tmp(6,:,:,:,:)
Qtot_flux		= Q_flux + Qpen_flux

Qtot_flux_diff		= Qtot_flux(1,:,:,:)-Qtot_flux(0,:,:,:)
copy_VarCoords(Q_flux,Qtot_flux)
copy_VarCoords_1(Q_flux(time|:,lat|:,lon|:,case|:),Qtot_flux_diff)
printVarSummary(Qtot_flux_diff)
printMinMax(Qtot_flux_diff,True)

vars_Terms(0,:,:,:,:)		= Qtot_flux
vars_TermsDiff(0,:,:,:)		= Qtot_flux_diff

ncdf = addfile("test-cal-not-terms-vertdiff3-Qpenflux.nc" ,"r")

; read mixed layer depth
HMXL				= ncdf->HMXL
HMXL_diff			= HMXL(1,:,:,:)-HMXL(0,:,:,:) ; (case,time,lat,lon)
copy_VarCoords_1(HMXL(time|:,lat|:,lon|:,case|:),HMXL_diff)

deltaT_mld			= ncdf->deltaT_mld
printMinMax(deltaT_mld,True)

; dhdt
dhdt				= HMXL
dhdt				= 0.
do itime=1,ntime-1
dhdt(:,itime,:,:)		= HMXL(:,itime,:,:)-HMXL(:,itime-1,:,:)
end do

dhdt_diff			= dhdt(1,:,:,:)-dhdt(0,:,:,:)
copy_VarCoords_1(dhdt(time|:,lat|:,lon|:,case|:),dhdt_diff)

h_dhdt				= (-1.)*dhdt/where(HMXL.ne.0,HMXL,HMXL@_FillValue)*deltaT_mld
copy_VarCoords(HMXL,h_dhdt)

h_dhdt_diff			= h_dhdt(1,:,:,:)-h_dhdt(0,:,:,:)
copy_VarCoords_1(h_dhdt(time|:,lat|:,lon|:,case|:),h_dhdt_diff)

vars_Terms(dimsizes(varTerms),:,:,:,:)			= h_dhdt
vars_TermsDiff(dimsizes(varTerms),:,:,:)		= h_dhdt_diff

do iterm=0,3 ; represent sum of all, no h_dhdt, no Qpenflux, no vertdiff3
; calculate the average of terms: "netFlux","u_dTdx","v_dTdy","hordiff","w_dTdz","vertdiff3","Qpenflux","h_dhdt"
vars_Terms_sum				= dim_sum_n_Wrap(vars_Terms(:dimsizes(varTerms)-iterm,:,:,:,:),0)
vars_TermsDiff_sum			= dim_sum_n_Wrap(vars_TermsDiff(:dimsizes(varTerms)-iterm,:,:,:),0)

vars_Terms(dimsizes(varTerms)+1+iterm,:,:,:,:)			= vars_Terms_sum
vars_TermsDiff(dimsizes(varTerms)+1+iterm,:,:,:)		= vars_TermsDiff_sum
end do

; read mixed layer mean Temperature data
SST				= ncdf->Temp_mldmean
SST_diff			= SST(1,:,:,:)-SST(0,:,:,:) ; (case,time,lat,lon)
copy_VarCoords_1(SST(time|:,lat|:,lon|:,case|:),SST_diff)

; dSSTdt
dSSTdt				= SST
dSSTdt				= 0.
do itime=1,ntime-1
dSSTdt(:,itime,:,:)		= SST(:,itime,:,:)-SST(:,itime-1,:,:)
end do

dSSTdt_diff			= dSSTdt(1,:,:,:)-dSSTdt(0,:,:,:)
copy_VarCoords_1(dSSTdt(time|:,lat|:,lon|:,case|:),dSSTdt_diff)

; residuals
residual			= dSSTdt
residual			= dSSTdt-vars_Terms(dimsizes(varTerms)+1,:,:,:,:)
printVarSummary(residual)

residual_diff			= residual(1,:,:,:)-residual(0,:,:,:)
copy_VarCoords_1(residual(time|:,lat|:,lon|:,case|:),residual_diff)

; copy dSSTdt and residual into vars_Terms arrays
vars_Terms(dimsizes(varTerms)+5,:,:,:,:)		= dSSTdt
vars_Terms(dimsizes(varTerms)+6,:,:,:,:)		= residual

vars_TermsDiff(dimsizes(varTerms)+5,:,:,:)		= dSSTdt_diff
vars_TermsDiff(dimsizes(varTerms)+6,:,:,:)		= residual_diff

printVarSummary(vars_Terms)
printVarSummary(vars_TermsDiff)
printMinMax(vars_TermsDiff(dimsizes(varTerms)-1,:,:,:),True)

print("nvars= "+nvars+" ncase= "+ncase+" nlat= "+nlat+" nlon= "+nlon)

; time integral
xtime			= fspan(0,ntime-1,ntime)

yval			= vars_Terms(vars|:,case|:,lat|:,lon|:,time|:)
yval_diff		= vars_TermsDiff(vars|:,lat|:,lon|:,time|:)

print(yval&lat)
print(yval&lon)

simpne1_tmp			= new((/nvars,ncase,nlat,nlon,ntime+1/),typeof(vars_Terms))
simpne1_tmp_Diff		= new((/nvars,nlat,nlon,ntime+1/),typeof(vars_Terms))

simpne1_tmp			= 0
simpne1_tmp_Diff		= 0

do itime=1,ntime
; attention: function "simpne" is only used when var numbers larger greater than 3.
;so I prescribe the first and second integral value manually.
	; for direct SST difference evolution, it is not time integral, I think.
	if(itime.eq.1)then
		simpne1_tmp(:,:,:,:,itime)		= yval(:,:,:,:,itime-1)
		simpne1_tmp_Diff(:,:,:,itime)		= yval_diff(:,:,:,itime-1)
	print("itime= "+itime+"  "+"simpne1= "+simpne1_tmp_Diff(dimsizes(varTerms)-1,{0},{180},itime))
		continue
	end if

	if(itime.eq.2)then
		simpne1_tmp(:,:,:,:,itime)		= simpne1_tmp(:,:,:,:,itime-1)+yval(:,:,:,:,itime-1)
		simpne1_tmp_Diff(:,:,:,itime)		= simpne1_tmp_Diff(:,:,:,itime-1)+yval_diff(:,:,:,itime-1)
	print("itime= "+itime+"  "+"simpne1= "+simpne1_tmp_Diff(dimsizes(varTerms)-1,{0},{180},itime))
		continue
	end if

	simpne1_tmp(:,:,:,:,itime)			= simpne(xtime(:itime-1),yval(:,:,:,:,:itime-1))
	simpne1_tmp_Diff(:,:,:,itime)			= simpne(xtime(:itime-1),yval_diff(:,:,:,:itime-1))
;	print("itime= "+itime+"  "+"simpne1= "+simpne1_tmp_Diff(dimsizes(varTerms)-1,{0},{180},itime))
	print("itime= "+itime+"  "+"simpne1= "+simpne1_tmp_Diff(5,{-16.5:-16.5},{178.5:179.5},itime))


end do

simpne1_tmp&time	= fspan(0,ntime,ntime+1)
simpne1_tmp_Diff&time	= fspan(0,ntime,ntime+1)

printVarSummary(simpne1_tmp)

simpne1			= simpne1_tmp(vars|:,case|:,time|:,lat|:,lon|:)
simpne1_Diff		= simpne1_tmp_Diff(vars|:,time|:,lat|:,lon|:)

; extremely strange values appear near 180E, 15S-20S, so I change them into missing values. 
; it's important to longitude-value plot. pay attention!

simpne1(:,:,:,{-16.5:-16.5},{178.5:179.5}) 		= -9999
simpne1_Diff(:,:,{-16.5:-16.5},{178.5:179.5}) 		= -9999

simpne1@_FillValue	= -9999
simpne1_Diff@_FillValue	= -9999

printVarSummary(simpne1)

system("/bin/rm -r simpne1-vertdiff3-Qpenflux-"+varTerms(dimsizes(varTerms)-1)+"-nca-island.nc")
ncdf = addfile("simpne1-vertdiff3-Qpenflux-"+varTerms(dimsizes(varTerms)-1)+"-nca-island.nc","c")
ncdf->simpne1 = simpne1
ncdf->simpne1_Diff = simpne1_Diff

end if ; DoTimeInteg


;####################################################
; read new file
;####################################################
ncin 			= addfile("simpne1-vertdiff3-Qpenflux-"+varTerms(dimsizes(varTerms)-1)+"-nca-island.nc","r")
simpne1 		= ncin->simpne1
simpne1_Diff 	 	= ncin->simpne1_Diff
printVarSummary(simpne1)
print(simpne1(0,0,:,{-10:-20},{175:180}))

dims			= dimsizes(simpne1)
nvars			= dims(0)
ncases			= dims(1)
ntime			= dims(2)
nlat			= dims(3)
nlon			= dims(4)
print(dims)


delete(dims)
;-----------------------------------------------
; Do region averages...
;-----------------------------------------------
;latSs				= (/-10,-10,-10,-20,-20,-20,-20,-20,-20,-5,-5,-5/)
;latEs				= (/-5,-5,-5,-10,-10,-10,-5,-5,-5,5,5,5/)
;lonSs				= (/180,230,260,180,230,260,180,230,260,240,200,160/)
;lonEs				= (/225,255,285,225,255,285,225,255,285,275,235,195/)

latSs	= (/-20,-20,-10,-20/)
latEs	= (/-5,-10,-5,-10/)
lonSs	= (/180,180,180,230/)
lonEs	= (/230,230,230,255/)

latS	= -20.
latE	= -15.

; calculate zonal mean over specific zonal bands
lat_here                        = simpne1&lat
rad                             = 4.0*atan(1.0)/180.
clat                            = cos(lat_here*rad)
clat!0                          = "lat"
clat&lat                        = lat_here

printVarSummary(simpne1)
printVarSummary(simpne1_Diff)

picktime			= 48

simpne1_zm_tmp			= dim_avg_wgt_n_Wrap(simpne1(:,:,picktime,{latS:latE},:),clat({latS:latE}),1,2)
simpne1_Diff_zm_tmp		= dim_avg_wgt_n_Wrap(simpne1_Diff(:,picktime,{latS:latE},:),clat({latS:latE}),1,1)

printVarSummary(simpne1_zm_tmp)
printVarSummary(simpne1_Diff_zm_tmp)

simpne1_zm			= simpne1_zm_tmp
simpne1_Diff_zm			= simpne1_Diff_zm_tmp

;==================================================================
; select needed varibles to be plotted
;==================================================================

simpne1_zm_need		=	simpne1_zm((/vector/),:,:)
simpne1_Diff_zm_need	=	simpne1_Diff_zm((/vector/),:)
printVarSummary(simpne1_zm_need)

; variables to be emphasized
simpne1_zm_empha		=	simpne1_zm((/vector_empha/),:,:)
simpne1_Diff_zm_empha	=	simpne1_Diff_zm((/vector_empha/),:)
printVarSummary(simpne1_zm_need)

legends_need			= 	legends_all((/vector/))
colors_need			=	colors_all((/vector/))
thick_need			=	thick_all((/vector/))
dashes_need			=	dashes_all((/vector/))

colors_empha			=	colors_all((/vector_empha/))
thick_empha			=	thick_all((/vector_empha/))
dashes_empha			=	dashes_all((/vector_empha/))


;==============================================
; Plotting starts here....
;==============================================

wks		 		= gsn_open_wks("x11","./figure/test-timeinteg-all-boxes-focused")

plot				= new(12,graphic)

res				= True

res@gsnFrame			= False
res@gsnDraw			= False

;res@vpWidthF			= 0.5
;res@vpHeightF			= 0.6

res@vpWidthF			= 0.7
res@vpHeightF			= 0.4


res@tiYAxisString		= ""
res@tiXAxisString		= ""

res@gsnYRefLine			= (/0.0/)

res@tmXBLabelFontHeightF        = 0.025
res@tmYLLabelFontHeightF        = 0.025
res@gsnStringFontHeightF	= 0.025

res@tiXAxisFontHeightF		= 0.030
res@tiYAxisFontHeightF		= 0.030



res1L			= res
res1L@xyLineColors	= colors_need
res1L@xyLineThicknesses	= thick_need
res1L@xyDashPatterns	= dashes_need


tmp			= (/abs(min(simpne1_zm_need)),abs(max(simpne1_zm_need))/)
res1L@trYMinF		= max(tmp)*(-1.)
res1L@trYMaxF		= max(tmp)
delete(tmp)

res1L@gsnLeftString	= "Park-RH"
plot(0)			= gsn_csm_xy(wks,simpne1_zm_need&lon,simpne1_zm_need(:,0,:),res1L)
res1L@gsnLeftString	= "Gauss-PDF"
plot(1)			= gsn_csm_xy(wks,simpne1_zm_need&lon,simpne1_zm_need(:,1,:),res1L)

; emphasize the relative values of dSST/dt, sumof4 and residual
res11			= res
res11@xyLineColors	= colors_empha
res11@xyLineThicknessF	= thick_empha
res11@xyDashPatterns	= dashes_empha

tmp			= (/abs(min(simpne1_zm_empha(:,0,:))),\
			    abs(max(simpne1_zm_empha(:,0,:)))/)
res11@trYMinF		= max(tmp)*(-1.)
res11@trYMaxF		= max(tmp)
delete(tmp)

plot(2)			= gsn_csm_xy(wks,simpne1_zm_need&lon,simpne1_zm_empha(:,0,:),res11)

; ------------------------
; Plot differences
; ------------------------
res2L			= res
res2L@xyLineColors	= colors_need
res2L@xyLineThicknesses	= thick_need
res2L@xyDashPatterns	= dashes_need

tmp			= (/abs(min(simpne1_Diff_zm_need)),abs(max(simpne1_Diff_zm_need))/)
res2L@trYMinF		= max(tmp)*(-1.)
res2L@trYMaxF		= max(tmp)
delete(tmp)


res2L@gsnLeftString	= "PDF-RH"
plot(3)			= gsn_csm_xy(wks,simpne1_zm_need&lon,simpne1_Diff_zm_need(:,:),res2L)

; emphasize the relative values of dSST/dt, sumof4 and residual
res22			= res
res22@xyLineColors	= colors_empha
res22@xyLineThicknessF	= thick_empha
res22@xyDashPatterns	= dashes_empha

tmp			= (/abs(min(simpne1_Diff_zm_empha(:,:))),\
			    abs(max(simpne1_Diff_zm_empha(:,:)))/)
res22@trYMinF		= max(tmp)*(-1.)
res22@trYMaxF		= max(tmp)
delete(tmp)

plot(4)			= gsn_csm_xy(wks,simpne1_zm_need&lon,simpne1_Diff_zm_empha(:,:),res22)

; add legends
legends_1		= legends_all((/vector1/))
legends_2		= legends_all((/vector2/))
legends_3		= legends_all

dashes_1		= dashes_all((/vector1/))
dashes_2		= dashes_all((/vector2/))

colors_1		= colors_all((/vector1/))
colors_2		= colors_all((/vector2/))

lgres                           = True
lgres@lgLabelFontHeightF        = 0.06
lgres@vpWidthF                  = 0.2
lgres@vpHeightF                 = 0.2
lgres@lgPerimOn                 = False

lgres1				= lgres
lgres1@lgDashIndexes		= dashes_1
lgres1@lgLineThicknessF		= 3.0
lgres1@lgLineColors		= colors_1
lgres1@lgItemOrder		= fspan(dimsizes(legends_1)-1,0,dimsizes(legends_1))

lgres2				= lgres
lgres2@lgDashIndexes		= dashes_2
lgres2@lgLineThicknessF		= 3.0
lgres2@lgLineColors		= colors_2
lgres2@lgItemOrder		= fspan(dimsizes(legends_2)-1,0,dimsizes(legends_2))

lgres3				= lgres
lgres3@vpWidthF                  = 0.2
lgres3@vpHeightF                 = 0.3
lgres3@lgDashIndexes		= dashes_all
;lgres3@lgLineThicknesses 	= thick_all ; looks like this one is problematic.
lgres3@lgLineThicknessF		= 3.0
lgres3@lgLineColors		= colors_all
lgres3@lgItemOrder		= fspan(dimsizes(legends_3)-1,0,dimsizes(legends_3))

lbid1                           = gsn_create_legend(wks,dimsizes(legends_1),legends_1,lgres1)
lbid2                           = gsn_create_legend(wks,dimsizes(legends_2),legends_2,lgres2)
lbid3                           = gsn_create_legend(wks,dimsizes(legends_3),legends_3,lgres3)

amres1				= True
amres1@amJust			= "TopLeft"
amres1@amParallelPosF		= -0.5
amres1@amOrthogonalPosF		= -0.5

amres2				= True
amres2@amJust			= "BottomLeft"
amres2@amParallelPosF		= -0.5
amres2@amOrthogonalPosF		= 0.5

amres3				= True
amres3@amJust			= "TopLeft"
amres3@amParallelPosF		= -0.5
amres3@amOrthogonalPosF		= -0.5

;if(iarea.eq.0)then
;gsn_legend_ndc(wks,dimsizes(legends_3),legends_3,0.1,0.95,lgres3)
;frame(wks)
;end if

;plot(5)				= gsn_csm_blank_plot(wks,res)
annoid1                         = gsn_add_annotation(plot(0),lbid1,amres1)
annoid2                         = gsn_add_annotation(plot(0),lbid2,amres2)
;annoid3                         = gsn_add_annotation(plot(3),lbid3,amres3)

;gsn_panel(wks,(/plot(5)/),(/1,1/),False)

;gsn_panel(wks,(/plot(2),plot(4)/),(/2,1/),False)
;gsn_panel(wks,(/plot(0),plot(1),plot(3)/),(/3,1/),False)

;gsn_panel(wks,(/plot(2),plot(4)/),(/1,2/),False)

resP				= True
resP@gsnPanelMainString		= ""
resP@amJust			= "BottomRight"
;if(iarea.eq.0)then
resP@gsnPanelFigureStrings	= (/"a)","b)","c)"/)
;else
;resP@gsnPanelFigureStrings	= (/"d)","e)","f)"/)
;end if

;gsn_panel(wks,(/plot(0),plot(1),plot(3)/),(/1,3/),resP)
;gsn_panel(wks,plot,(/2,3/),resP)
draw(plot(3))
frame(wks)

exit

;print(simpne1)




end do; do is=
                               
end
