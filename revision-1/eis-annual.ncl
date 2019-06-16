
function cal_eis(temp[*][*][*]:numeric,qv[*][*][*]:numeric,lev[*]:numeric)
begin

; Start to calculate EIS...
; 1. define constant parameters
RV 			= 461.5 ; J/kg/K
RD 			= 287.04 ; J/kg/K
EPS 			= RD/RV
CPD 			= 1005.7 ; J/kg/K
CPV 			= 1870.0 ; J/kg/K
P0 			= 1000.0 ; hPa
LV 			= 2501000 ; J/kg
gravit 			= 9.8 ; m/s2

; 1. calculate potential temperature ------------------------------------
P 			= conform(temp,lev,0)	; hPa
theta			= temp*(P0/P)^0.286 ; K
copy_VarCoords(temp,theta)
printVarSummary(theta)

; 2. calculate moist-adiabatic potential temp. gradient Gamma_m ---------------------
; 2.1 convert spc humidty to mixing ratio first
mr			= qv/(1.-qv)
copy_VarCoords(qv,mr)
printVarSummary(mr)

; 2.2 get saturation mixing ratio
; water vapor pressure
ev			 = mr*P/(EPS+mr) 
copy_VarCoords(mr,ev)
printVarSummary(ev)
; saturation water vapor pressure (hPa)
es			=6.112*exp(17.67*(temp-273.15)/(243.5+(temp-273.15)))
copy_VarCoords(temp,es)
printVarSummary(es)
; saturation mixing ratio
mrs			= EPS*es/(P-es)
copy_VarCoords(es,mrs)
printVarSummary(mrs)
; saturation spc humidity
qs			= mrs/(1.+mrs)
copy_VarCoords(mrs,qs)
printVarSummary(qs)

; 2.3 get moist-adiabatic potential temperature gradient 
Gamma		= gravit/CPD*(1.0-(1.0+LV*qs/RD/temp)/(1.0+LV^2*qs/CPD/RV/temp^2))
copy_VarCoords(qs,Gamma)
printVarSummary(Gamma)

; 3. calculate lifting condensation level: Bolton (1980) Eq.(21) ----------------------------
TC			= 2840.0/(3.5*log(temp(29,:,:))-log(ev(29,:,:))-4.805)+55.0
copy_VarCoords(temp(29,:,:),TC)

plcl		= P(29,:,:)*((TC/temp(29,:,:))^3.48)
copy_VarCoords(temp(29,:,:),plcl)

printMinMax(TC,True)
printMinMax(plcl,False)

; 4. get z700 and zlcl ---------------------------------------------------------------------
z700		= (RD*temp(29,:,:)/gravit)*log(P0/700.)
zlcl 		= (RD*temp(29,:,:)/gravit)*log(P0/plcl)

copy_VarCoords(temp(29,:,:),z700)
copy_VarCoords(temp(29,:,:),zlcl)

; 5. calculate estimated inversion strength (EIS) ------------------------------------------
LTS			= theta({700},:,:)-theta(29,:,:)
Gamma_850		= Gamma({850},:,:)
dz			= z700-zlcl
EIS			= where(dz.lt.0,LTS,LTS-Gamma({850},:,:)*dz)
copy_VarCoords(temp(0,:,:),LTS)
copy_VarCoords(temp(0,:,:),EIS)

return ([/LTS,EIS/])

end


begin 

fdir				= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/atm/"

case 				= (/"BC5_f19g16_cosp",\
		 		    "BC5_f19g16_mac2_cosp"/)

append				= "_ANN_climo.nc"

vars				= (/"T","Q"/)

nlat				= 96
nlon				= 144
nlev			 	= 30


f1				= addfile(fdir+case(0)+append,"r")

lev				= f1->lev

var				= new((/dimsizes(case),dimsizes(vars),nlev,nlat,nlon/),float)
var!0				= "cases"
var&cases			= fspan(0,dimsizes(case)-1,dimsizes(case))
var!1				= "vars"
var&vars			= fspan(0,dimsizes(vars)-1,dimsizes(vars))

var				= var

do icase=0,dimsizes(case)-1
	ftmp1			= addfile(fdir+case(icase)+append,"r")
	do ivar=0,dimsizes(vars)-1
		var(icase,ivar,:,:,:)	 = ftmp1->$vars(ivar)$(0,:,:,:)
	end do ; do ivar=
	lev			= ftmp1->lev
end do ; do icase=

printVarSummary(lev)


LTS				= new((/dimsizes(case),nlat,nlon/),double)
LTS!0				= "cases"
LTS&cases			= fspan(0,dimsizes(case)-1,dimsizes(case))

EIS				= LTS

do icase=0,dimsizes(case)-1

result			= cal_eis(var(icase,0,:,:,:),var(icase,1,:,:,:),lev)
LTS(icase,:,:)		= result[0]
EIS(icase,:,:)		= result[1]
end do ; do icase=

printVarSummary(LTS)
printMinMax(LTS,True)
printMinMax(EIS,True)


; get the difference
LTS_diff			= LTS(1,:,:)-LTS(0,:,:)
copy_VarCoords(LTS(0,:,:),LTS_diff)
EIS_diff			= EIS(1,:,:)-EIS(0,:,:)
copy_VarCoords(EIS(0,:,:),EIS_diff)

;LTS_diff			= where(LTS_diff.gt.0,LTS_diff,LTS@_FillValue)
;EIS_diff			= where(EIS_diff.gt.0,EIS_diff,LTS@_FillValue)
printVarSummary(EIS_diff)
printMinMax(EIS_diff,True)

; -----------------------------------------------------------------------
;create plot
; -----------------------------------------------------------------------
varout				= (/"EIS",	"LTS"/)
vars_units			= (/"K",	"K"/)
EIS_levs			= (/-3.,-2,-1,0,1,2,3,4,5,6,8,10,12/)
LTS_levs			= (/5.,6,7,8,9,10,12,14,16,18,20,22,24/)

EIS_dlevs_1			= (/-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1/)
EIS_dlevs_2			= (/-2.0,-1.6,-1.2,-0.8,-0.4,-0.2,0,0.2,0.4,0.8,1.2,1.6,2.0/)

LTS_dlevs_1			= (/-1.,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1/)
LTS_dlevs_2			= (/-2.0,-1.6,-1.2,-0.8,-0.4,-0.2,0,0.2,0.4,0.8,1.2,1.6,2.0/)


vars_levs			= (/EIS_levs,		LTS_levs/)
vars_dlevs_1			= (/EIS_dlevs_1,	LTS_dlevs_1/)
vars_dlevs_2			= (/EIS_dlevs_2,	LTS_dlevs_2/)

varout_data			= new((/dimsizes(varout),dimsizes(case),nlat,nlon/),double)
varout_data_diff		= new((/dimsizes(varout),nlat,nlon/),double)

varout_data!0			= "vars"
varout_data&vars		= fspan(0,dimsizes(varout)-1,dimsizes(varout))
varout_data(0,:,:,:)		= EIS
varout_data(1,:,:,:)		= LTS

varout_data_diff!0		= "vars"
varout_data_diff&vars		= fspan(0,dimsizes(varout)-1,dimsizes(varout))
varout_data_diff(0,:,:)		= EIS_diff
varout_data_diff(1,:,:)		= LTS_diff

printVarSummary(varout_data)

do ivar=0,dimsizes(varout)-1

wks 				= gsn_open_wks("pdf","./figure/"+varout(ivar))

cmap1 					= read_colormap_file("WhiteBlueGreenYellowRed")
;cmap2 					= read_colormap_file("BlueRed")
cmap2					= read_colormap_file("BlueDarkRed18")
;cmap2					= read_colormap_file("nrl_sirkes")
;cmap2					= read_colormap_file("hotcold_18lev")
;cmap2					= read_colormap_file("BlueWhiteOrangeRed")


res                      		= True
res@gsnDraw 				= False
res@gsnFrame 				= False
;res@mpProjection			= "Hammer"
;res@mpPerimOn              		= False	 

res@cnFillDrawOrder             	= "Predraw"

res@cnFillOn            		= True
res@cnLinesOn           		= False
res@cnLineLabelsOn      		= False
res@cnFillPalette       		= cmap1(:240,:)

res@mpCenterLonF         		= 180.

fontheight				= 0.025
res@gsnStringFontHeightF		= fontheight
res@tmXBLabelFontHeightF		= fontheight
res@tmYLLabelFontHeightF		= fontheight

;res@tmXBOn			 	= False
;res@tmYLOn				= False
;res@tmXTOn				= False
;res@tmYROn				= False
res@tmXBTickSpacingF			= 60.
res@tmYLTickSpacingF			= 45.

;res@gsnCenterString			= varout(ivar)
res@gsnRightString			= vars_units(ivar)

;res@lbOrientation			= "vertical"
res@lbLabelFontHeightF			= fontheight
res@pmLabelBarOrthogonalPosF 		= .25           ; move whole thing down

res@cnLevelSelectionMode        	= "ExplicitLevels"
res@cnLevels                    	= vars_levs(ivar,:)

plot					= new(3,graphic)

res@gsnLeftString			= "DEF"
plot(0)					= gsn_csm_contour_map(wks,varout_data(ivar,0,:,:),res)

res@gsnLeftString			= "NEW"
plot(1)					= gsn_csm_contour_map(wks,varout_data(ivar,1,:,:),res)

delete(res@cnFillPalette)
res@cnLevels                    	= vars_dlevs_2(ivar,:)
res@cnFillPalette       		= cmap2
res@gsnLeftString			= "NEW-DEF"
plot(2)					= gsn_csm_contour_map(wks,varout_data_diff(ivar,:,:),res)

resP = True
resP@gsnPanelMainString = varout(ivar)

resP@gsnPanelXWhiteSpacePercent		= 2.5
resP@gsnPanelYWhiteSpacePercent		= 5.

gsn_panel(wks,plot,(/2,2/),resP)

delete(res)

end do ; do ivar=

end
