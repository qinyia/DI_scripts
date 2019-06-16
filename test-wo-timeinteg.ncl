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

;varTerms			= (/"netFlux","u_dTdx","v_dTdy","entrain"/)
;varTerms			= (/"netFlux","u_dTdx","v_dTdy","w_dTdz","entrain","hordiff","vertdiff3","Qpenflux"/)
varTerms			= (/"netFlux","u_dTdx","v_dTdy","hordiff","w_dTdz","entrain","vertdiff3","Qpenflux"/)

;varTerms			= (/"netFlux","u_dTdx","v_dTdy","entrain","vertdiff3"/)
;varTerms			= (/"netFlux","u_dTdx","v_dTdy","entrain","vertdiff2","vertdiff3"/)
;varTerms			= (/"netFlux","u_dTdx","v_dTdy","entrain","vertdiff0","vertdiff2"/)
;varTerms			= (/"netFlux","UET_mldmean","VNT_mldmean","WTT_mldmean","vertdiff"/)
;varTerms			= (/"netFlux","u_dTdx","v_dTdy","entrain","hordiff"/)

;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff"/)
varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","hordiff_diff","w_dTdz_diff","entrain_diff","vertdiff3_diff","Qpenflux_diff"/)
;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff","vertdiff3_diff"/)
;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff","vertdiff2_diff","vertdiff3_diff"/)
;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff","vertdiff0_diff","vertdiff2_diff"/)
;varTermsDiff			= (/"netFlux_diff","UET_mldmean_diff","VNT_mldmean_diff","WTT_mldmean_diff","vertdiff_diff"/)
;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff","hordiff_diff"/)

colors				= (/"brown4","aquamarine","cadetblue","blue","pink2","magenta2","purple2","brown1"/)
;dashes				= (/       0,      0,      0,      0,      0,        0,           0,       0/)
dashes				= fspan(0,0,dimsizes(varTerms))
;thick				= (/       3,      3,      3,      3,      3,        3,           3,       3/)
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
ncdf = addfile("test-cal-terms-vertdiff3-Qpenflux.nc" ,"r")

dSizes                          = getfiledimsizes(ncdf)  

; netFlux,u_dTdx,v_dTdy,entrain,sum,dSSTdt,residual
vars_Terms         	= new((/dimsizes(varTerms)+3,dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_Terms!0		= "vars"
vars_Terms&vars		= ispan(0,dimsizes(varTerms)-1+3,1)

vars_TermsDiff         	= new((/dimsizes(varTermsDiff)+3,dSizes(1),dSizes(2),dSizes(3)/),float)
vars_TermsDiff!0	= "vars"
vars_TermsDiff&vars	= ispan(0,dimsizes(varTermsDiff)-1+3,1)

do ivar=0,dimsizes(varTerms)-1
        vars_Terms(ivar,:,:,:,:)           =ncdf->$varTerms(ivar)$ 
if(varTerms(ivar).eq."vertdiff0")then
	vars_Terms(ivar,:,:,:,:)	= vars_Terms(ivar,:,:,:,:)*(-1.)
end if

end do

do ivar=0,dimsizes(varTermsDiff)-1
	vars_TermsDiff(ivar,:,:,:)		= ncdf->$varTermsDiff(ivar)$
if(varTerms(ivar).eq."vertdiff0_diff")then
	vars_TermsDiff(ivar,:,:,:)	= vars_TermsDiff(ivar,:,:,:)*(-1.)
end if


end do

nvars			= dimsizes(vars_Terms&vars)
ncase			= dimsizes(vars_Terms&case)
nlat			= dimsizes(vars_Terms&lat)
nlon			= dimsizes(vars_Terms&lon)
ntime			= dimsizes(vars_Terms&time)

; calculate the average of four terms: netFlux, u_dTdx, v_dTdy, entrain
vars_Terms_sum			= dim_sum_n_Wrap(vars_Terms(:dimsizes(varTerms)-1,:,:,:,:),0)
vars_TermsDiff_sum		= dim_sum_n_Wrap(vars_TermsDiff(:dimsizes(varTerms)-1,:,:,:),0)

vars_Terms(dimsizes(varTerms),:,:,:,:)			= vars_Terms_sum
vars_TermsDiff(dimsizes(varTerms),:,:,:)		= vars_TermsDiff_sum

; read mixed layer mean Temperature data
ncdf = addfile("test-cal-not-terms-vertdiff3-Qpenflux.nc" ,"r")

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
residual			= dSSTdt-vars_Terms_sum
printVarSummary(residual)
printMinMax(residual,True)

residual_diff			= residual(1,:,:,:)-residual(0,:,:,:)
copy_VarCoords_1(residual(time|:,lat|:,lon|:,case|:),residual_diff)
printMinMax(residual_diff,True)


; copy dSSTdt and residual into vars_Terms arrays

vars_Terms(dimsizes(varTerms)+1,:,:,:,:)		= dSSTdt
vars_Terms(dimsizes(varTerms)+2,:,:,:,:)		= residual

vars_TermsDiff(dimsizes(varTerms)+1,:,:,:)		= dSSTdt_diff
vars_TermsDiff(dimsizes(varTerms)+2,:,:,:)		= residual_diff

printVarSummary(vars_Terms)
printVarSummary(vars_TermsDiff)
printMinMax(vars_TermsDiff(dimsizes(varTerms)-1,:,:,:),True)

print("nvars= "+nvars+" ncase= "+ncase+" nlat= "+nlat+" nlon= "+nlon)


dims			= dimsizes(vars_Terms)
nvars			= dims(0)
ncases			= dims(1)
ntime			= dims(2)
nlat			= dims(3)
nlon			= dims(4)
delete(dims)

;plot tests

;plot			=new(nvars,graphic)
;terms			= array_append_record (varTerms, (/"sumof4","dSSTdt","residual"/), 0)
;
;
;wks		 	= gsn_open_wks("x11","ocn-heat-budget-evo-each-term-check-1220-1")
;res			= True
;
;res@gsnFrame		= False
;res@gsnDraw		= False
;
;res@vpWidthF		= 0.75
;res@vpHeightF		= 0.4
;
;res@cnFillOn		= True
;res@cnLinesOn		= False
;res@cnLineLabelsOn	= False
;
;res@gsnLeftString	= ""
;
;res@cnLevelSelectionMode        = "ManualLevels"
;res@cnMinLevelValF              = -5.
;res@cnMaxLevelValF              = 5.
;res@cnLevelSpacingF             = 1. 
;
;
;do itime=1,ntime-1
;	res@gsnCenterString		= itime+1
;do iplot=0,dimsizes(plot)-1
;	res@gsnRightString		= terms(iplot)
;	plot(iplot)			= gsn_csm_contour(wks,vars_Terms(iplot,0,itime,:,:),res)
;end do
;gsn_panel(wks,plot,(/3,3/),False)
;end do



;-----------------------------------------------
; Do region averages...
;-----------------------------------------------
latSs				= (/-10,-10,-10,-20,-20,-20,-20,-20,-20,-5,-5,-5/)
latEs				= (/-5,-5,-5,-10,-10,-10,-5,-5,-5,5,5,5/)
lonSs				= (/180,230,260,180,230,260,180,230,260,240,200,160/)
lonEs				= (/225,255,285,225,255,285,225,255,285,275,235,195/)

region				= new((/dimsizes(latSs)/),float)
region				= fspan(0,dimsizes(latSs)-1,dimsizes(latSs))



lat_here                        = vars_Terms&lat
rad                             = 4.0*atan(1.0)/180.
clat                            = cos(lat_here*rad)
clat!0                          = "lat"
clat&lat                        = lat_here

wks		 	= gsn_open_wks("pdf","./figure/test-wo-timeinteg-all-boxes")

res			= True

res@gsnFrame		= False
res@gsnDraw		= False

res@vpWidthF		= 0.75
res@vpHeightF		= 0.4

res@tiYAxisString	= ""
res@tiXAxisString	= ""

res@gsnYRefLine		= (/0.0/)


do iarea = 0,dimsizes(latSs)-1

simpne1_areamean_tmp		= wgt_areaave_Wrap(vars_Terms(:,:,:,{latSs(iarea):latEs(iarea)},{lonSs(iarea):lonEs(iarea)}),clat({latSs(iarea):latEs(iarea)}),1.0,0)

simpne1_Diff_areamean_tmp	= wgt_areaave_Wrap(vars_TermsDiff(:,:,{latSs(iarea):latEs(iarea)},{lonSs(iarea):lonEs(iarea)}),clat({latSs(iarea):latEs(iarea)}),1.0,0)

runnum				= 12

simpne1_areamean		=	runave_n_Wrap(simpne1_areamean_tmp,runnum,0,2)
simpne1_Diff_areamean		=	runave_n_Wrap(simpne1_Diff_areamean_tmp,runnum,0,1)

printVarSummary(simpne1_areamean)
printVarSummary(simpne1_Diff_areamean)
printMinMax(simpne1_Diff_areamean(dimsizes(varTerms)-2,:),True)


plot			= new(12,graphic)

;nsteps			= 256/nvars
;colors			= ispan(2,256,nsteps)
colors_all		=  array_append_record(colors,(/"grey","black","grey"/),0)
dashes_all		=  array_append_record(dashes,(/0.,0.,14./),0)
thick_all		=  array_append_record(thick,(/3.,3.,3./),0)

res1L			= res
res1L@xyLineColors	= colors_all
res1L@xyLineThicknesses	= thick_all
res1L@xyDashPatterns	= dashes_all


plot(0)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_areamean(:,0,:),res1L)
plot(1)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_areamean(:,1,:),res1L)

res2L			= res
res2L@xyLineColors	= colors_all
res2L@xyLineThicknesses	= thick_all
res2L@xyDashPatterns	= dashes_all


res2L@gsnRightString	= "latS="+latSs(iarea)+" latE="+latEs(iarea)+" lonS="+lonSs(iarea)+" lonE="+lonEs(iarea)


plot(2)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_Diff_areamean(:,:),res2L)

; emphasize the relative values of dSST/dt, sumof4 and residual
res22			= res
res22@xyLineColors	= colors_all(dimsizes(varTerms):)
res22@xyLineThicknessF	= thick_all(dimsizes(varTerms):)
res22@xyDashPatterns	= dashes_all(dimsizes(varTerms):)


plot(3)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_Diff_areamean(dimsizes(varTerms):,:),res22)

; add legends
;legends_1				= array_append_record (varTerms, (/"sumof4","dSSTdt","residual"/), 0)
;legends_2				= array_append_record (varTerms, (/"sumof4","dSSTdt","residual"/), 0)
legends_1				= varTerms(:5)
legends_2				= array_append_record(varTerms(6:),(/"sumof4","dSSTdt","residual"/),0)

;(/"netFlux","u_dTdx","v_dTdy","entrain","hordiff","vertdiff"/)

lgres                           = True
lgres@lgLabelFontHeightF        = 0.06
lgres@vpWidthF                  = 0.20
lgres@vpHeightF                 = 0.15
lgres@lgPerimOn                 = False

lgres1				= lgres
lgres1@lgDashIndexes		= dashes_all(:5)
lgres1@lgLineThicknessF		= thick_all(:5)
lgres1@lgLineColors		= colors_all(:5)

lgres2				= lgres
lgres2@lgDashIndexes		= dashes_all(6:)
lgres2@lgLineThicknessF 	= thick_all(6:)
lgres2@lgLineColors		= colors_all(6:)


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


annoid1                         = gsn_add_annotation(plot(0),lbid1,amres1)
annoid2                         = gsn_add_annotation(plot(0),lbid2,amres2)

draw(plot(0))
frame(wks)


;gsn_panel(wks,plot(0:1),(/2,1/),False)
gsn_panel(wks,plot(0:3),(/2,2/),False)

end do 
exit

;print(simpne1)





end do; do is=
                               
end
