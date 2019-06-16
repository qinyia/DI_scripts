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

varTerms			= (/"HMXL","deltaT_mld","VDC_T_mldmean","dTdz_mld"/)
varTermsDiff			= (/"HMXL_diff","deltaT_mld","VDC_T_mldmean","dTdz_mld"/)
colors				= (/"blue","red3","green2","orange"/)

;varTerms			= (/"HMXL","VDC_T_mldmean","dTdz_mld"/)
;varTermsDiff			= (/"HMXL_diff","VDC_T_mldmean","dTdz_mld"/)
;colors				= (/"blue","green2","orange"/)



dashes				= fspan(0,0,dimsizes(varTerms)+1)
thick				= fspan(4,4,dimsizes(varTerms)+1)

compare_type 			= (/"cmip-obs"/)
case				= (/"def","new"/)

;Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"
Dir			= "/home/lyl/WORK4/cesm1_2_1/archive/"

compset			= (/"ocn"/)

addnums			= 2+2		; dhdt,h_dhdt,kdTdz,kdTdz_h

do is = 0,dimsizes(season)-1

;========================================================
; read NC file
;========================================================
ncdf = addfile("test-cal-not-terms-vertdiff3-Qpenflux.nc" ,"r")

dSizes                          = getfiledimsizes(ncdf)  
print(dSizes)

vars_Terms         	= new((/dimsizes(varTerms)+addnums,dSizes(0),dSizes(1),dSizes(2),dSizes(3)/),float)
vars_Terms!0		= "vars"
vars_Terms&vars		= ispan(0,dimsizes(varTerms)-1+addnums,1)

vars_TermsDiff         	= new((/dimsizes(varTermsDiff)+addnums,dSizes(1),dSizes(2),dSizes(3)/),float)
vars_TermsDiff!0	= "vars"
vars_TermsDiff&vars	= ispan(0,dimsizes(varTermsDiff)-1+addnums,1)

;--------------------------------------------------------
; read mixed layer depth
;--------------------------------------------------------
do ivar=0,dimsizes(varTerms)-1
        vars_Terms(ivar,:,:,:,:)           =ncdf->$varTerms(ivar)$ 
end do

vars_TermsDiff(:,:,:,:)			= vars_Terms(:,1,:,:,:)-vars_Terms(:,0,:,:,:)
copy_VarCoords_1(vars_Terms(vars|:,time|:,lat|:,lon|:,case|:),vars_TermsDiff)
printVarSummary(vars_TermsDiff)

nvars			= dimsizes(vars_Terms&vars)
ncase			= dimsizes(vars_Terms&case)
nlat			= dimsizes(vars_Terms&lat)
nlon			= dimsizes(vars_Terms&lon)
ntime			= dimsizes(vars_Terms&time)


HMXL						= vars_Terms(0,:,:,:,:)
HMXL_diff					= HMXL(1,:,:,:)-HMXL(0,:,:,:)
copy_VarCoords_1(HMXL(time|:,lat|:,lon|:,case|:),HMXL_diff)
printVarSummary(vars_Terms)

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

h_dhdt				= (-1.)*dhdt/where(HMXL.ge.0,HMXL,HMXL@_FillValue)
copy_VarCoords(HMXL,h_dhdt)
printVarSummary(h_dhdt)
printMinMax(h_dhdt,True)

h_dhdt_diff			= h_dhdt(1,:,:,:)-h_dhdt(0,:,:,:)
copy_VarCoords_1(h_dhdt(time|:,lat|:,lon|:,case|:),h_dhdt_diff)
printVarSummary(h_dhdt_diff)
printMinMax(h_dhdt_diff,True)

;VDC_T_mld*dTdz_mld

VDC_T_mld			= vars_Terms(2,:,:,:,:)
dTdz_mld			= vars_Terms(3,:,:,:,:)

kdTdz				= VDC_T_mld*dTdz_mld*(-1.)
copy_VarCoords(VDC_T_mld,kdTdz)
printVarSummary(kdTdz)

kdTdz_diff			= kdTdz(1,:,:,:)-kdTdz(0,:,:,:)
copy_VarCoords_1(kdTdz(time|:,lat|:,lon|:,case|:),kdTdz_diff)
printVarSummary(kdTdz_diff)

kdTdz_h				= kdTdz/where(HMXL.ge.0,HMXL,HMXL@_FillValue)
copy_VarCoords(VDC_T_mld,kdTdz_h)

kdTdz_h_diff			= kdTdz_h(1,:,:,:)-kdTdz_h(0,:,:,:)
copy_VarCoords_1(kdTdz_h(time|:,lat|:,lon|:,case|:),kdTdz_h_diff)
printVarSummary(kdTdz_h_diff)

vars_Terms(dimsizes(varTerms):dimsizes(varTerms)+addnums-1,:,:,:,:)		= (/dhdt,h_dhdt,kdTdz,kdTdz_h/)
vars_TermsDiff(dimsizes(varTerms):dimsizes(varTerms)+addnums-1,:,:,:)		= (/dhdt_diff,h_dhdt_diff,kdTdz_diff,kdTdz_h_diff/)


nvars			= dimsizes(vars_Terms&vars)
ncase			= dimsizes(vars_Terms&case)
nlat			= dimsizes(vars_Terms&lat)
nlon			= dimsizes(vars_Terms&lon)
ntime			= dimsizes(vars_Terms&time)

varTerms_1		= array_append_record(varTerms,(/"dhdt","h_dhdt","kdTdz","kdTdz_h"/),0)


; test spatial pattern of mixed-layer depth
;wks			= gsn_open_wks("x11","test")
;res			= True
;
;res@cnFillOn		= True
;res@cnLinesOn		= False
;res@cnLineLabelsOn	= False
;
;res@cnLevelSelectionMode        = "ManualLevels"
;res@cnMinLevelValF              = -50.
;res@cnMaxLevelValF              = 50.
;res@cnLevelSpacingF             = 5.
;
;do itime=0,12
;plot	= gsn_csm_contour(wks,HMXL_diff(itime,:,:),res)
;end do


;-----------------------------------------------
; Do region averages...
;-----------------------------------------------
;latSs				= (/-10,-10,-10,-20,-20,-20,-20,-20,-20,-5,-5,-5/)
;latEs				= (/-5,-5,-5,-10,-10,-10,-5,-5,-5,5,5,5/)
;lonSs				= (/180,230,260,180,230,260,180,230,260,240,200,160/)
;lonEs				= (/225,255,285,225,255,285,225,255,285,275,235,195/)

latSs	= (/-20,-10/)
latEs	= (/-10,-5/)
lonSs	= (/180,180/)
lonEs	= (/230,230/)

lat_here                        = vars_Terms&lat
rad                             = 4.0*atan(1.0)/180.
clat                            = cos(lat_here*rad)
clat!0                          = "lat"
clat&lat                        = lat_here

runnum				= 12

wks		 	= gsn_open_wks("pdf","./figure/test-HMXL-1-runum-"+runnum)

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


simpne1_areamean		=	runave_n_Wrap(simpne1_areamean_tmp,runnum,0,2)
simpne1_Diff_areamean		=	runave_n_Wrap(simpne1_Diff_areamean_tmp,runnum,0,1)

printVarSummary(simpne1_areamean)
printVarSummary(simpne1_Diff_areamean)
printMinMax(simpne1_Diff_areamean(dimsizes(varTerms)-1,:),True)


plot			= new((/nvars,3/),graphic)

do ivar = 0,nvars-1

res@gsnLeftString	= varTerms_1(ivar)


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

tmp			= (/abs(min(simpne1_areamean(ivar,:,:))),abs(max(simpne1_areamean(ivar,:,:)))/)
res1L@trYMinF		= max(tmp)*(-1.)
res1L@trYMaxF		= max(tmp)
;res1L@trYMinF		= min(simpne1_areamean(ivar,:,:))
;res1L@trYMaxF		= max(simpne1_areamean(ivar,:,:))

delete(tmp)

plot(ivar,0)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_areamean(ivar,0,:),res1L)
plot(ivar,1)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_areamean(ivar,1,:),res1L)

; ------------------------
; Plot differences
; ------------------------
res2L			= res
res2L@xyLineColors	= colors_all
res2L@xyLineThicknesses	= thick_all
res2L@xyDashPatterns	= dashes_all

tmp			= (/abs(min(simpne1_Diff_areamean(ivar,:))),abs(max(simpne1_Diff_areamean(ivar,:)))/)
res2L@trYMinF		= max(tmp)*(-1.)
res2L@trYMaxF		= max(tmp)
;res2L@trYMinF		= min(simpne1_Diff_areamean(ivar,:))
;res2L@trYMaxF		= max(simpne1_Diff_areamean(ivar,:))
delete(tmp)

res2L@gsnRightString	= "latS="+latSs(iarea)+" latE="+latEs(iarea)+" lonS="+lonSs(iarea)+" lonE="+lonEs(iarea)

plot(ivar,2)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_Diff_areamean(ivar,:),res2L)

if(False)then
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

;if(iarea.eq.0)then
;gsn_legend_ndc(wks,dimsizes(legends_3),legends_3,0.1,0.95,lgres3)
;frame(wks)
;
;plot(5)				= gsn_csm_blank_plot(wks,res)
;annoid1                         = gsn_add_annotation(plot(5),lbid1,amres1)
;annoid2                         = gsn_add_annotation(plot(5),lbid2,amres2)
;
;gsn_panel(wks,(/plot(5)/),(/1,1/),False)
;end if

end if

;gsn_panel(wks,(/plot(ivar,:)/),(/3,1/),False)


end do ; do ivar=

dimsplot		= dimsizes(plot)
gsn_panel(wks,ndtooned(plot),(/dimsplot(0),3/),False)

end do  ; do iarea=
exit

;print(simpne1)




end do; do is=
                               
end
