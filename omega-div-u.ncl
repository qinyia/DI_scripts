load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

function sfvp(lev[*]:numeric,lat[*]:numeric,lon[*]:numeric,\
              U_1[*][*][*]:numeric,V_1[*][*][*]:numeric,OMEGA_1[*][*][*],\
	      gw[*]:numeric,latS:float,latE:float)

begin

nlev 			= dimsizes(lev)
nlat 			= dimsizes(lat)
nlon 			= dimsizes(lon)

sfvp_1 			= new((/2,nlev,nlat,nlon/),float)

sfvp_1(0,:,:,:) 	= U_1
sfvp_1 			= uv2sfvpG(U_1,V_1)

sf_1 			= sfvp_1(0,:,:,:)		; stream function -- non-divergent
vp_1 			= sfvp_1(1,:,:,:)		; velocity potential -- irrotation


DvpDx_1 		= new((/nlev,nlat,nlon/),float)
DvpDx_1 		= vp_1
;DvpDx_1 		= doubletofloat(center_finite_diff_n(vp_1,lon*111000,False,0,2))

DvpDx_1 		= doubletofloat(center_finite_diff_n(vp_1,lat*111000,False,0,1))

printVarSummary(DvpDx_1)

DsfDy_1 		= new((/nlev,nlat,nlon/),float)
DsfDy_1 		= sf_1
DsfDy_1 		= -1.*doubletofloat(center_finite_diff_n(sf_1,lon*111000,False,0,2))
printVarSummary(DsfDy_1)


DvpDx_1_mm 		= dim_avg_wgt_n_Wrap(DvpDx_1(:,{latS:latE},:),gw({latS:latE}),0,1)
lev111 			= DvpDx_1_mm&lev
lev111@units 		= "hPa"
DvpDx_1_mm&lev 		= lev111 
printVarSummary(DvpDx_1_mm)

DsfDy_1_mm 		= dim_avg_wgt_n_Wrap(DsfDy_1(:,{latS:latE},:),gw({latS:latE}),0,1)
lev111 			= DsfDy_1_mm&lev
lev111@units 		= "hPa"
DsfDy_1_mm&lev 		= lev111 
printVarSummary(DsfDy_1_mm)


OMEGA_1_mm 		= dim_avg_wgt_n_Wrap(OMEGA_1(:,{latS:latE},:),gw({latS:latE}),0,1)
lev111 			= OMEGA_1_mm&lev
lev111@units 		= "hPa"
OMEGA_1_mm&lev 		= lev111 
printVarSummary(OMEGA_1_mm)

printMinMax(OMEGA_1_mm,False)

OMEGA_1_mm 		= (/-1.0*OMEGA_1_mm*1e2/)

return([/DvpDx_1_mm,DsfDy_1_mm,OMEGA_1_mm/])

end


begin

;***************file 1********************
;fdir 			= "/home/lyl/WORK3/qinyi/double-ITCZ/data/amip/"
;case 			= (/"FC5_f19f19_MG15_cosp","FC5_f19f19_MG15_mac3_cosp"/)

season			= (/"ANN","SON"/)

compare_type 		= (/"amip-obs","cmip-obs"/)
;compare_type 		= (/"cmip-obs"/)

Dir			= "/home/lyl/WORK3/qinyi/double-ITCZ/data/"

do icomp = 0,dimsizes(compare_type)-1
print(compare_type(icomp))

if(compare_type(icomp) .eq."cmip-obs")then
	case	 	= (/"BC5_f19g16","BC5_f19g16_mac2"/)
	fdir 		= (/Dir+"atm/",Dir+"atm/"/)
end if
if(compare_type(icomp) .eq."amip-obs")then
	case		= (/"FC5_f19f19_MG15_cosp","FC5_f19f19_MG15_mac3_cosp"/)
	fdir 		= (/Dir+"amip/",Dir+"amip/"/)
end if

do is			= 0, dimsizes(season)-1

append 			= "_"+season(is)+"_climo.nc"

f0 			= addfile(fdir(0)+case(0)+append,"r")
f1 			= addfile(fdir(1)+case(1)+append,"r")

gw 			= f1->gw
lev 			= f1->lev
lat 			= f1->lat
lon 			= f1->lon

nlat 			= dimsizes(lat)
nlon 			= dimsizes(lon)
nlev 			= dimsizes(lev)


U_1 			= f0->U(0,:,:,:)
V_1 			= f0->V(0,:,:,:)
OMEGA_1			= f0->OMEGA(0,:,:,:)

U_2 			= f1->U(0,:,:,:)
V_2 			= f1->V(0,:,:,:)
OMEGA_2 		= f1->OMEGA(0,:,:,:)


;1111111111111111111111111111111

latS 			= -10.
latE 			= -5.

lonS			= 120.
lonE			= 280.

result 			= sfvp(lev,lat,lon,U_1,V_1,OMEGA_1,gw,latS,latE)
DvpDx_1_mm 		= result[0]
DsfDy_1_mm 		= result[1]
OMEGA_1_mm 		= result[2]
delete(result)

printVarSummary(DvpDx_1_mm)

result 			= sfvp(lev,lat,lon,U_2,V_2,OMEGA_2,gw,latS,latE)
DvpDx_2_mm 		= result[0]
DsfDy_2_mm 		= result[1]
OMEGA_2_mm 		= result[2]
delete(result)

;************************create plot***************
wks = gsn_open_wks("pdf","./figure/lon-height-omega-div-irrot-u_"+compare_type(icomp)+"_"+season(is)+"_"+latS+"to"+latE)
gsn_define_colormap(wks,"BlueRed")



vcres 				= True
vcres@gsnDraw 			= False
vcres@gsnFrame 			= False

vcres@trYReverse 		= True

vcres@trXMinF 			= lonS
vcres@trXMaxF 			= lonE

vcres@trYMinF 			= 100.

vcres@vpHeightF 		= 0.3
vcres@vpWidthF 			= 0.7

vcres@tiYAxisString 		= "Pressure (hPa)"
vcres@gsnCenterString		= season(is)

vcres@cnLinesOn 		= False
vcres@cnLineLabelsOn 		= False
vcres@cnFillOn 			= True
vcres@cnInfoLabelOn 		= False

maxval                          = max(DvpDx_1_mm)
minval                          = min(DvpDx_1_mm)

tmp                             = (/abs(maxval),abs(minval)/)
exmval                          = max(tmp)
data_levels                     = fspan(-1.*exmval,exmval,11)

vcres@cnLevelSelectionMode	= "ExplicitLevels"
vcres@cnLevels           	= data_levels

vcres@cnLineLabelFormat           = "0@;*.3f"                      ; Only one value after decimal point
vcres@lbLabelStrings              = sprintf("%5.1f",data_levels)   ; Format the labelbar labels

vcres@vcRefMagnitudeF  		= 3.0
vcres@vcRefLengthF     		= 0.045                     ; what the ref length is
vcres@vcGlyphStyle     		= "CurlyVector" 
vcres@vcLineArrowThicknessF 	= 3.5
vcres@vcRefAnnoOrthogonalPosF 	= -1.15
vcres@vcRefAnnoArrowSpaceF 	= 2.0
vcres@vcRefAnnoString2On 	= False
vcres@vcVectorDrawOrder		= "Postdraw"

vcres@vcMinDistanceF		= 0.040

ovplot1				= new(3, graphic)
ovplot2 			= new(3, graphic)

DvpDx_21_mm 			= DvpDx_1_mm 
DvpDx_21_mm 			= (/DvpDx_2_mm-DvpDx_1_mm/)

DsfDy_21_mm 			= DsfDy_1_mm 
DsfDy_21_mm 			= (/DsfDy_2_mm-DsfDy_1_mm/)

OMEGA_21_mm 			= OMEGA_1_mm 
OMEGA_21_mm 			= (/OMEGA_2_mm-OMEGA_1_mm/)


vcres@gsnLeftString = "DEF"
vcres@gsnRightString 		= "W and divergent U"
ovplot1(0) = gsn_csm_pres_hgt_vector(wks,DvpDx_1_mm, DvpDx_1_mm,OMEGA_1_mm,vcres)
vcres@gsnRightString 		= "W and rotational U"
ovplot2(0) = gsn_csm_pres_hgt_vector(wks,DsfDy_1_mm, DsfDy_1_mm,OMEGA_1_mm,vcres)

vcres@gsnLeftString = "NEW"
vcres@gsnRightString 		= "W and divergent U"
ovplot1(1) = gsn_csm_pres_hgt_vector(wks,DvpDx_2_mm, DvpDx_2_mm,OMEGA_2_mm,vcres)
vcres@gsnRightString 		= "W and rotational U"
ovplot2(1) = gsn_csm_pres_hgt_vector(wks,DsfDy_2_mm, DsfDy_2_mm,OMEGA_2_mm,vcres)

if(compare_type(icomp).eq."amip-obs")then
	vcres@vcRefMagnitudeF  = 1.
end if

if(compare_type(icomp).eq."cmip-obs")then
	vcres@vcRefMagnitudeF  = 2.
end if
vcres@gsnLeftString = "NEW-DEF"
vcres@gsnRightString 		= "W and divergent U"
ovplot1(2) = gsn_csm_pres_hgt_vector(wks,DvpDx_21_mm, DvpDx_21_mm,OMEGA_21_mm,vcres)
vcres@gsnRightString 		= "W and rotational U"
ovplot2(2) = gsn_csm_pres_hgt_vector(wks,DsfDy_21_mm, DsfDy_21_mm,OMEGA_21_mm,vcres)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, ovplot1, (/3,1/), pres)
gsn_panel(wks, ovplot2, (/3,1/), pres)

ovplot	= (/ovplot1,ovplot2/)
gsn_panel(wks,ndtooned(ovplot),(/2,3/),pres)

end do ; do is=

end do ; do icomp=
end

