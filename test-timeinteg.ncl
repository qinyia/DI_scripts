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
varTerms			= (/"netFlux","u_dTdx","v_dTdy","w_dTdz","entrain","hordiff","vertdiff3","Qpenflux"/)
;varTerms			= (/"netFlux","u_dTdx","v_dTdy","entrain","vertdiff3"/)
;varTerms			= (/"netFlux","u_dTdx","v_dTdy","entrain","vertdiff2","vertdiff3"/)
;varTerms			= (/"netFlux","u_dTdx","v_dTdy","entrain","vertdiff0","vertdiff2"/)
;varTerms			= (/"netFlux","UET_mldmean","VNT_mldmean","WTT_mldmean","vertdiff"/)
;varTerms			= (/"netFlux","u_dTdx","v_dTdy","entrain","hordiff"/)

;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff"/)
varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","w_dTdz_diff","entrain_diff","hordiff_diff","vertdiff3_diff","Qpenflux_diff"/)
;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff","vertdiff3_diff"/)
;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff","vertdiff2_diff","vertdiff3_diff"/)
;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff","vertdiff0_diff","vertdiff2_diff"/)
;varTermsDiff			= (/"netFlux_diff","UET_mldmean_diff","VNT_mldmean_diff","WTT_mldmean_diff","vertdiff_diff"/)
;varTermsDiff			= (/"netFlux_diff","u_dTdx_diff","v_dTdy_diff","entrain_diff","hordiff_diff"/)

compare_type 			= (/"cmip-obs"/)
case				= (/"def","new"/)

;DoTimeInteg			= True
DoTimeInteg			= False


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
end do

do ivar=0,dimsizes(varTermsDiff)-1
	vars_TermsDiff(ivar,:,:,:)		= ncdf->$varTermsDiff(ivar)$
end do

nvars			= dimsizes(vars_Terms&vars)
ncase			= dimsizes(vars_Terms&case)
nlat			= dimsizes(vars_Terms&lat)
nlon			= dimsizes(vars_Terms&lon)
ntime			= dimsizes(vars_Terms&time)

; calculate the sum of four terms: netFlux, u_dTdx, v_dTdy, entrain
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


;-----------------------------------------------
; Do time integrals...
;-----------------------------------------------
if(DoTimeInteg)then

xtime			= fspan(0,ntime-1,ntime)

yval			= vars_Terms(vars|:,case|:,lat|:,lon|:,time|:)
yval_diff		= vars_TermsDiff(vars|:,lat|:,lon|:,time|:)

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
	print("itime= "+itime+"  "+"simpne1= "+simpne1_tmp_Diff(dimsizes(varTerms)-1,{0},{180},itime))

end do

simpne1_tmp&time	= fspan(0,ntime,ntime+1)
simpne1_tmp_Diff&time	= fspan(0,ntime,ntime+1)

printVarSummary(simpne1_tmp)

simpne1			= simpne1_tmp(vars|:,case|:,time|:,lat|:,lon|:)
simpne1_Diff		= simpne1_tmp_Diff(vars|:,time|:,lat|:,lon|:)

simpne1@_FillValue	= -9999
simpne1_Diff@_FillValue	= -9999

printVarSummary(simpne1)
;printMinMax(simpne1_Diff(5,:,:,:),True)

system("/bin/rm -r simpne1-vertdiff3-Qpenflux.nc")
ncdf = addfile("simpne1-vertdiff3-Qpenflux.nc","c")
ncdf->simpne1 = simpne1
ncdf->simpne1_Diff = simpne1_Diff

end if ; DoTimeInteg


ncin 			= addfile("simpne1-vertdiff3-Qpenflux.nc","r")
simpne1 		= ncin->simpne1
simpne1_Diff 	 	= ncin->simpne1_Diff

dims			= dimsizes(simpne1)
nvars			= dims(0)
ncases			= dims(1)
ntime			= dims(2)
nlat			= dims(3)
nlon			= dims(4)
delete(dims)

;-----------------------------------------------
; Do region averages...
;-----------------------------------------------
;latSs				= (/-15,-20,-15,-5,-5,-5/)
;latEs				= (/-5,-5,-5,5,5,5/)
;lonSs				= (/180,250,210,240,200,160/)
;lonEs				= (/220,280,250,270,230,190/)

;latSs				= (/-15,-20,-15,-5,-5,-5/)
;latEs				= (/-10,-10,-18,5,5,5/)
;lonSs				= (/180,250,210,240,200,160/)
;lonEs				= (/220,280,250,270,230,190/)

latSs				= (/-10,-10,-10,-20,-20,-20,-20,-20,-20,-5,-5,-5/)
latEs				= (/-5,-5,-5,-10,-10,-10,-5,-5,-5,5,5,5/)
lonSs				= (/180,230,260,180,230,260,180,230,260,240,200,160/)
lonEs				= (/225,255,285,225,255,285,225,255,285,275,235,195/)



;latSs				= (/-20/)
;latEs				= (/-5/)
;lonSs				= (/180/)
;lonEs				= (/230/)
region				= new((/dimsizes(latSs)/),float)
region				= fspan(0,dimsizes(latSs)-1,dimsizes(latSs))



lat_here                                = simpne1&lat
rad                                     = 4.0*atan(1.0)/180.
clat                                    = cos(lat_here*rad)
clat!0                                  = "lat"
clat&lat                                = lat_here


wks		 	= gsn_open_wks("x11","./figure/test-timeinteg-all-boxes")

res			= True

res@gsnFrame		= False
res@gsnDraw		= False

res@vpWidthF		= 0.75
res@vpHeightF		= 0.4

res@tiYAxisString	= ""
res@tiXAxisString	= "time(month)"


do iarea = 0,dimsizes(latSs)-1

simpne1_areamean_tmp		= wgt_areaave_Wrap(simpne1(:,:,:,{latSs(iarea):latEs(iarea)},{lonSs(iarea):lonEs(iarea)}),clat({latSs(iarea):latEs(iarea)}),1.0,0)

simpne1_Diff_areamean_tmp	= wgt_areaave_Wrap(simpne1_Diff(:,:,{latSs(iarea):latEs(iarea)},{lonSs(iarea):lonEs(iarea)}),clat({latSs(iarea):latEs(iarea)}),1.0,0)

runnum				= 3
simpne1_areamean		=	runave_n_Wrap(simpne1_areamean_tmp,runnum,0,2)
simpne1_Diff_areamean		=	runave_n_Wrap(simpne1_Diff_areamean_tmp,runnum,0,1)

printVarSummary(simpne1_areamean)
printVarSummary(simpne1_Diff_areamean)
printMinMax(simpne1_Diff_areamean(dimsizes(varTerms)-1,:),True)



plot			= new(4,graphic)

nsteps			= 256/nvars
colors			= ispan(2,256,nsteps)

res1L			= res
res1L@xyLineColors	= colors
res1L@xyLineThicknessF	= 3.0
res1L@xyDashPatterns	= fspan(0,0,nvars)


plot(0)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_areamean(:,1,:),res1L)


res2L			= res
res2L@xyLineColors	= colors
res2L@xyLineThicknessF	= 3.0
res2L@xyDashPatterns	= fspan(0,0,nvars)

;res2L@trYMinF		= -20.
;res2L@trYMaxF		= 20.
res2L@gsnYRefLine	= (/0.0/)

res2L@gsnRightString	= "latS="+latSs(iarea)+" latE="+latEs(iarea)+" lonS="+lonSs(iarea)+" lonE="+lonEs(iarea)


plot(1)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_Diff_areamean(:,:),res2L)

res22			= res
res22@xyLineColors	= colors(dimsizes(varTerms):)
res22@xyLineThicknessF	= 3.0
res22@xyDashPatterns	= fspan(0,0,nvars)

res22@gsnYRefLine	= (/0.0/)


plot(2)			= gsn_csm_xy(wks,simpne1_areamean&time,simpne1_Diff_areamean(dimsizes(varTerms):,:),res22)


; add legends
legends_1				= array_append_record (varTerms, (/"sumof4","dSSTdt","residual"/), 0)
legends_2				= array_append_record (varTerms, (/"sumof4","dSSTdt","residual"/), 0)

;(/"netFlux","u_dTdx","v_dTdy","entrain","hordiff","vertdiff"/)

lgres                           = True
lgres@lgLabelFontHeightF        = 0.06
lgres@vpWidthF                  = 0.20
lgres@vpHeightF                 = 0.20
lgres@lgPerimOn                 = False

lgres1				= lgres
lgres1@lgDashIndexes		= res1L@xyDashPatterns
lgres1@lgLineThicknessF 	= res1L@xyLineThicknessF
lgres1@lgLineColors		= res1L@xyLineColors

lgres2				= lgres
lgres2@lgDashIndexes		= res2L@xyDashPatterns
lgres2@lgLineThicknessF 	= res2L@xyLineThicknessF
lgres2@lgLineColors		= res2L@xyLineColors


lbid1                           = gsn_create_legend(wks,dimsizes(legends_1),legends_1,lgres1)
lbid2                           = gsn_create_legend(wks,dimsizes(legends_2),legends_2,lgres2)

amres                           = True
amres@amJust			= "TopLeft"
amres@amParallelPosF            = -0.50
amres@amOrthogonalPosF          = -0.50

;annoid1                         = gsn_add_annotation(plot(0),lbid1,amres)
annoid2                         = gsn_add_annotation(plot(2),lbid2,amres)

;gsn_panel(wks,plot(0:1),(/2,1/),False)
gsn_panel(wks,plot(0:2),(/2,2/),False)

end do 
exit

;print(simpne1)



end do; do is=
                               
end
