load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

load "./functions_lat_interpolate.ncl"
load "./functions_get_resource.ncl"
load "./functions_get_var.ncl"

begin

;*********************************************
; read in data
;*********************************************

;season=(/"ANN","DJF","MAM","JJA","SON"/)

season				=(/"ANN"/)

varTimeDepthLatLon		=(/"TEMP","UVEL","VVEL","WVEL"/)
varTimeLatLon			=(/"HMXL","SHF"/)

varTerms			= (/"netFlux","u_dTdx","v_dTdy","hordiff","w_dTdz","vertdiff3","Qpenflux","shflx","lhflx","lwnet","shf_qsw"/)
;					0        1        2         3        4          5          6         7       8        9      10
varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","hordiff_diff","w_dTdz_diff","vertdiff3_diff","Qpenflux_diff","shflx_diff","lhflx_diff","lwnet_diff","shf_qsw_diff"/)
colors				= (/"brown4","aquamarine","cadetblue","blue","magenta2","purple2","brown1","seagreen1","orange","yellow","cyan"/)


dashes				= fspan(0,0,dimsizes(varTerms))
thick				= fspan(2,2,dimsizes(varTerms))

compare_type 			= (/"cmip-obs"/)
case				= (/"def","new"/)

;Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"
Dir			= "/home/lyl/WORK4/cesm1_2_1/archive/"

compset			= (/"ocn"/)


do is = 0,dimsizes(season)-1

;========================================================
; read NC file
;========================================================
ncdf = addfile("test-cal-terms-vertdiff3-Qpenflux-R1.nc" ,"r")

dSizes                          = getfiledimsizes(ncdf)  

; netFlux,u_dTdx,v_dTdy,entrain,sum,dSSTdt,residual
vars_Terms         	= new((/dimsizes(varTerms)+4+3,dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_Terms!0		= "vars"
vars_Terms&vars		= ispan(0,dimsizes(varTerms)-1+4+3,1)

vars_TermsDiff         	= new((/dimsizes(varTermsDiff)+4+3,dSizes(1),dSizes(2),dSizes(3)/),float)
vars_TermsDiff!0	= "vars"
vars_TermsDiff&vars	= ispan(0,dimsizes(varTermsDiff)-1+4+3,1)

do ivar=0,dimsizes(varTerms)-1
        vars_Terms(ivar,:,:,:,:)           =ncdf->$varTerms(ivar)$ 
end do

do ivar=0,dimsizes(varTermsDiff)-1
	vars_TermsDiff(ivar,:,:,:)		= ncdf->$varTermsDiff(ivar)$
end do

nvars			= dimsizes(vars_Terms&vars)
ncase			= dimsizes(vars_Terms&case)
nlat			= dimsizes(vars_Terms&lat)
nlon			= dimsizes(vars_Terms&lon)
ntime			= dimsizes(vars_Terms&time)


ncdf = addfile("test-cal-not-terms-vertdiff3-Qpenflux-R1.nc" ,"r")
;--------------------------------------------------------
; read mixed layer depth
;--------------------------------------------------------
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
printVarSummary(dhdt)
printMinMax(dhdt,True)

dhdt_diff			= dhdt(1,:,:,:)-dhdt(0,:,:,:)
copy_VarCoords_1(dhdt(time|:,lat|:,lon|:,case|:),dhdt_diff)
printVarSummary(dhdt_diff)
printMinMax(dhdt,True)

h_dhdt				= (-1.)*dhdt/where(HMXL.ne.0,HMXL,HMXL@_FillValue)*deltaT_mld
copy_VarCoords(HMXL,h_dhdt)
printVarSummary(h_dhdt)
printMinMax(h_dhdt,True)

h_dhdt_diff			= h_dhdt(1,:,:,:)-h_dhdt(0,:,:,:)
copy_VarCoords_1(h_dhdt(time|:,lat|:,lon|:,case|:),h_dhdt_diff)
printVarSummary(h_dhdt_diff)
printMinMax(h_dhdt_diff,True)

vars_Terms(dimsizes(varTerms),:,:,:,:)			= h_dhdt
vars_TermsDiff(dimsizes(varTerms),:,:,:)		= h_dhdt_diff

do iterm=0,3 ; represent sum of all, no h_dhdt, no Qpenflux, no vertdiff3
; calculate the average of terms: "netFlux","u_dTdx","v_dTdy","hordiff","w_dTdz","vertdiff3","Qpenflux","h_dhdt"
vars_Terms_sum				= dim_sum_n_Wrap(vars_Terms(:dimsizes(varTerms)-iterm,:,:,:,:),0)
vars_TermsDiff_sum			= dim_sum_n_Wrap(vars_TermsDiff(:dimsizes(varTerms)-iterm,:,:,:),0)

vars_Terms(dimsizes(varTerms)+1+iterm,:,:,:,:)			= vars_Terms_sum
vars_TermsDiff(dimsizes(varTerms)+1+iterm,:,:,:)		= vars_TermsDiff_sum
end do
;-------------------------------------------
; read mixed layer mean Temperature data
;-------------------------------------------
SST				= ncdf->Temp_mldmean
SST_diff			= SST(1,:,:,:)-SST(0,:,:,:) ; (case,time,lat,lon)
copy_VarCoords_1(SST(time|:,lat|:,lon|:,case|:),SST_diff)

; dSSTdt
dSSTdt				= SST
dSSTdt				= 0.
do itime=1,ntime-1
dSSTdt(:,itime,:,:)		= SST(:,itime,:,:)-SST(:,itime-1,:,:)
end do
printVarSummary(dSSTdt)
printMinMax(dSSTdt,True)

dSSTdt_diff			= dSSTdt(1,:,:,:)-dSSTdt(0,:,:,:)
copy_VarCoords_1(dSSTdt(time|:,lat|:,lon|:,case|:),dSSTdt_diff)
printVarSummary(dSSTdt)
printMinMax(dSSTdt,True)

; residuals
residual			= dSSTdt
;residual			= dSSTdt-vars_Terms_sum
residual			= dSSTdt-vars_Terms(dimsizes(varTerms)+1,:,:,:,:)
printVarSummary(residual)
printMinMax(residual,True)

residual_diff			= residual(1,:,:,:)-residual(0,:,:,:)
copy_VarCoords_1(residual(time|:,lat|:,lon|:,case|:),residual_diff)
printMinMax(residual_diff,True)


; copy dSSTdt and residual into vars_Terms arrays

vars_Terms(dimsizes(varTerms)+5,:,:,:,:)		= dSSTdt
vars_Terms(dimsizes(varTerms)+6,:,:,:,:)		= residual

vars_TermsDiff(dimsizes(varTerms)+5,:,:,:)		= dSSTdt_diff
vars_TermsDiff(dimsizes(varTerms)+6,:,:,:)		= residual_diff

printVarSummary(vars_Terms)
printVarSummary(vars_TermsDiff)
printMinMax(vars_TermsDiff(dimsizes(varTerms)-1,:,:,:),True)

print("nvars= "+nvars+" ncase= "+ncase+" nlat= "+nlat+" nlon= "+nlon)
printMinMax(vars_TermsDiff(9,:,:,:),True)

dims			= dimsizes(vars_Terms)
nvars			= dims(0)
ncases			= dims(1)
ntime			= dims(2)
nlat			= dims(3)
nlon			= dims(4)
delete(dims)


;-----------------------------------------------
; Do region averages...
;-----------------------------------------------
;latSs				= (/-10,-10,-10,-20,-20,-20,-20,-20,-20,-5,-5,-5/)
;latEs				= (/-5,-5,-5,-10,-10,-10,-5,-5,-5,5,5,5/)
;lonSs				= (/180,230,260,180,230,260,180,230,260,240,200,160/)
;lonEs				= (/225,255,285,225,255,285,225,255,285,275,235,195/)

latSs	= (/-20,-20/)
latEs	= (/-5,-10/)
lonSs	= (/180,230/)
lonEs	= (/230,255/)


region				= new((/dimsizes(latSs)/),float)
region				= fspan(0,dimsizes(latSs)-1,dimsizes(latSs))



lat_here                        = vars_Terms&lat
rad                             = 4.0*atan(1.0)/180.
clat                            = cos(lat_here*rad)
clat!0                          = "lat"
clat&lat                        = lat_here

wks		 	= gsn_open_wks("pdf","./figure/test-wo-timeinteg-all-boxes-"+varTerms(dimsizes(varTerms)-1)+"-no-legends-R1")

res			= True

res@gsnFrame		= False
res@gsnDraw		= False

res@vpWidthF		= 0.75
res@vpHeightF		= 0.4

res@tiYAxisString	= ""
res@tiXAxisString	= ""

res@gsnYRefLine		= (/0.0/)

do iarea = 0,dimsizes(latSs)-1

res@gsnRightString	= "latS="+latSs(iarea)+" latE="+latEs(iarea)+" lonS="+lonSs(iarea)+" lonE="+lonEs(iarea)

simpne1_areamean_tmp		= wgt_areaave_Wrap(vars_Terms(:,:,:,{latSs(iarea):latEs(iarea)},{lonSs(iarea):lonEs(iarea)}),clat({latSs(iarea):latEs(iarea)}),1.0,0)

simpne1_Diff_areamean_tmp	= wgt_areaave_Wrap(vars_TermsDiff(:,:,{latSs(iarea):latEs(iarea)},{lonSs(iarea):lonEs(iarea)}),clat({latSs(iarea):latEs(iarea)}),1.0,0)

runnum				= 1

simpne1_areamean		=	runave_n_Wrap(simpne1_areamean_tmp,runnum,0,2)
simpne1_Diff_areamean		=	runave_n_Wrap(simpne1_Diff_areamean_tmp,runnum,0,1)

printVarSummary(simpne1_areamean)
printVarSummary(simpne1_Diff_areamean)

print(simpne1_Diff_areamean(9,:))
print(varTerms(9))

exit
printMinMax(simpne1_Diff_areamean(dimsizes(varTerms)-2,:),True)


plot			= new(12,graphic)

;nsteps			= 256/nvars
;colors			= ispan(2,256,nsteps)
colors_all		=  array_append_record(colors,(/"green",\
							"grey","burlywood3","bisque2","antiquewhite1",\
							"black","grey"/),0)
dashes_all		=  array_append_record(dashes,(/0.,\
							0.,0.,0.,0.,\
							0.,14./),0)
thick_all		=  array_append_record(thick,(/2.,\
						       4.,3.,3.,3.,\
						       4.,4./),0)

res1L			= res
res1L@xyLineColors	= colors_all
res1L@xyLineThicknesses	= thick_all
res1L@xyDashPatterns	= dashes_all

res1L@gsnRightString	= "latS="+latSs(iarea)+" latE="+latEs(iarea)+" lonS="+lonSs(iarea)+" lonE="+lonEs(iarea)

tmp			= (/abs(min(simpne1_areamean)),abs(max(simpne1_areamean))/)
res1L@trYMinF		= max(tmp)*(-1.)
res1L@trYMaxF		= max(tmp)
delete(tmp)

plot(0)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_areamean(:,0,:),res1L)
plot(1)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_areamean(:,1,:),res1L)

; emphasize the relative values of dSST/dt, sumof4 and residual
res11			= res
res11@xyLineColors	= colors_all(dimsizes(varTerms)+1:)
res11@xyLineThicknessF	= thick_all(dimsizes(varTerms)+1:)
res11@xyDashPatterns	= dashes_all(dimsizes(varTerms)+1:)

tmp			= (/abs(min(simpne1_areamean(dimsizes(varTerms)+1:,0,:))),\
			    abs(max(simpne1_areamean(dimsizes(varTerms)+1:,0,:)))/)
res11@trYMinF		= max(tmp)*(-1.)
res11@trYMaxF		= max(tmp)
delete(tmp)

plot(2)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_areamean(dimsizes(varTerms)+1:,0,:),res11)

; ------------------------
; Plot differences
; ------------------------
res2L			= res
res2L@xyLineColors	= colors_all
res2L@xyLineThicknesses	= thick_all
res2L@xyDashPatterns	= dashes_all

tmp			= (/abs(min(simpne1_Diff_areamean)),abs(max(simpne1_Diff_areamean))/)
res2L@trYMinF		= max(tmp)*(-1.)
res2L@trYMaxF		= max(tmp)
delete(tmp)

res2L@gsnRightString	= "latS="+latSs(iarea)+" latE="+latEs(iarea)+" lonS="+lonSs(iarea)+" lonE="+lonEs(iarea)

plot(3)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_Diff_areamean(:,:),res2L)

; emphasize the relative values of dSST/dt, sumof4 and residual
res22			= res
res22@xyLineColors	= colors_all(dimsizes(varTerms)+1:)
res22@xyLineThicknessF	= thick_all(dimsizes(varTerms)+1:)
res22@xyDashPatterns	= dashes_all(dimsizes(varTerms)+1:)

tmp			= (/abs(min(simpne1_Diff_areamean(dimsizes(varTerms)+1:,:))),abs(max(simpne1_Diff_areamean(dimsizes(varTerms)+1:,:)))/)
res22@trYMinF		= max(tmp)*(-1.)
res22@trYMaxF		= max(tmp)
delete(tmp)
plot(4)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_Diff_areamean(dimsizes(varTerms)+1:,:),res22)


; add legends
legends_1				= varTerms(:dimsizes(varTerms)-2)
legends_2				= array_append_record(varTerms(dimsizes(varTerms)-1:),\
							      (/"h_dhdt",\
							        "sumof4","no h_dhdt","no_Qpenflux","no_vertdiff",\
								"dSSTdt","residual"/),0)

legends_3				=  array_append_record(varTerms,\
							      (/"h_dhdt",\
							        "sumof4","no h_dhdt","no_Qpenflux","no_vertdiff",\
								"dSSTdt","residual"/),0)
;(/"netFlux","u_dTdx","v_dTdy","entrain","hordiff","vertdiff"/)

lgres                           = True
lgres@lgLabelFontHeightF        = 0.06
lgres@vpWidthF                  = 0.20
lgres@vpHeightF                 = 0.15
lgres@lgPerimOn                 = False

lgres1				= lgres
lgres1@lgDashIndexes		= dashes_all(:dimsizes(varTerms)-2)
;lgres1@lgLineThicknesses	= thick_all(:dimsizes(varTerms)-2)
lgres1@lgLineThicknessF		= 3.0
lgres1@lgLineColors		= colors_all(:dimsizes(varTerms)-2)

lgres2				= lgres
lgres2@lgDashIndexes		= dashes_all(dimsizes(varTerms)-1:)
;lgres2@lgLineThicknesses	= thick_all(dimsizes(varTerms)-1:)
lgres2@lgLineThicknessF		= 3.0
lgres2@lgLineColors		= colors_all(dimsizes(varTerms)-1:)

lgres3				= lgres
lgres3@vpWidthF                  = 0.4
lgres3@vpHeightF                 = 0.7
lgres3@lgDashIndexes		= dashes_all
;lgres3@lgLineThicknesses 	= thick_all
lgres3@lgLineThicknessF		= 3.0
lgres3@lgLineColors		= colors_all

lbid1                           = gsn_create_legend(wks,dimsizes(legends_1),legends_1,lgres1)
lbid2                           = gsn_create_legend(wks,dimsizes(legends_2),legends_2,lgres2)

amres1                           = True
amres1@amJust			= "BottomLeft"
amres1@amParallelPosF            = -0.5
amres1@amOrthogonalPosF          = 1.1

amres2                           = True
amres2@amJust			= "BottomRight"
amres2@amParallelPosF            = 0.5
amres2@amOrthogonalPosF          = 1.1

if(iarea.eq.0)then
gsn_legend_ndc(wks,dimsizes(legends_3),legends_3,0.1,0.95,lgres3)
frame(wks)

plot(5)				= gsn_csm_blank_plot(wks,res)
annoid1                         = gsn_add_annotation(plot(5),lbid1,amres1)
annoid2                         = gsn_add_annotation(plot(5),lbid2,amres2)

gsn_panel(wks,(/plot(5)/),(/1,1/),False)
end if

gsn_panel(wks,(/plot(2),plot(4)/),(/2,1/),False)
gsn_panel(wks,(/plot(0),plot(1),plot(3)/),(/3,1/),False)

end do 
exit

;print(simpne1)




end do; do is=
                               
end
