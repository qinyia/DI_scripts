
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

;========================================================================
;========================================================================
function get_var_atm_TimeLatLon_addfiles(Dir:string,fdir_append:string,case:string,var:string,type:string)
begin

fils 			= systemfunc ("ls "+Dir+case+fdir_append+"*.nc") ; file paths
f_tmp     		= addfiles (fils, "r")                         ; add file

fils_2			= systemfunc ("ls "+Dir+case+fdir_append+"*0001-01.nc") ; file paths
f_tmp_2			= addfile(fils_2,"r")

ListSetType (f_tmp, "cat")

if(type.eq."model")then
	if(isfilevar(f_tmp_2,var))then
		var_tmp    			= f_tmp[:]->$var$
		if(var.eq."TAUX".or.var.eq."TAUY")then
			var_tmp = var_tmp*(-1.)
		end if
		if(var.eq."PSL")then
			var_tmp	= var_tmp/100. ;change from Pa to hPa
		end if
	else
		if(var.eq."PRECT")then
			PRECC_tmp   = f_tmp[:]->PRECC
			PRECL_tmp   = f_tmp[:]->PRECL
			
			var_tmp = PRECC_tmp
			var_tmp@long_name = "total precipitation rate"
			var_tmp   = PRECC_tmp+PRECL_tmp
			
			var_tmp = var_tmp*8.64e7
			var_tmp@units = "mm/day"
			
			delete([/PRECC_tmp,PRECL_tmp/])
		end if
		if(var.eq."NetSfcHeat")then
			FSNS_tmp 	= f_tmp[:]->FSNS
			FLNS_tmp 	= f_tmp[:]->FLNS
			LHFLX_tmp 	= f_tmp[:]->LHFLX
			SHFLX_tmp 	= f_tmp[:]->SHFLX
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = FSNS_tmp-FLNS_tmp-LHFLX_tmp-SHFLX_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp,LHFLX_tmp,SHFLX_tmp/])
		end if
		if(var.eq."NetAtmHeat")then
			FSNS_tmp 	= f_tmp[:]->FSNS
			FLNS_tmp 	= f_tmp[:]->FLNS
			LHFLX_tmp 	= f_tmp[:]->LHFLX
			SHFLX_tmp 	= f_tmp[:]->SHFLX
			FSNT_tmp 	= f_tmp[:]->FSNT
			FLUT_tmp 	= f_tmp[:]->FLNT

			PRECL_tmp	= f_tmp[:]->PRECL
			PRECC_tmp	= f_tmp[:]->PRECC
			PRECSC_tmp	= f_tmp[:]->PRECSC
			PRECSL_tmp 	= f_tmp[:]->PRECSL
			QFLX_tmp	= f_tmp[:]->QFLX
			
			Lv = 2.501e6                   ; J/kg
			Lf = 3.337e5                   ; J/kg
			LHFLX_new = (/(Lv+Lf)*QFLX_tmp - Lf*1.e3*(PRECC_tmp+PRECL_tmp-PRECSC_tmp-PRECSL_tmp)/)
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = (FSNT_tmp-FLUT_tmp)-(FSNS_tmp-FLNS_tmp-LHFLX_new-SHFLX_tmp)
			var_tmp@units = "W/m2"
			delete([/FSNT_tmp,FLUT_tmp,FSNS_tmp,FLNS_tmp,LHFLX_tmp,SHFLX_tmp/])
		end if
		if(var.eq."NetToaHeat")then
			FSNT_tmp 	= f_tmp[:]->FSNT
			FLUT_tmp 	= f_tmp[:]->FLNT

			var_tmp = FSNT_tmp
			var_tmp@long_name = "Net TOA heat flux"
			var_tmp = (/FSNT_tmp-FLUT_tmp/)
			var_tmp@units = "W/m2"
			delete([/FSNT_tmp,FLUT_tmp/])
		end if
		if(var.eq."CrfToa")then
			FSNT_tmp	= f_tmp[:]->FSNT
			FSNTC_tmp	= f_tmp[:]->FSNTC
			FLNT_tmp	= f_tmp[:]->FLNT
			FLNTC_tmp	= f_tmp[:]->FLNTC

			var_tmp		= FSNT_tmp
			var_tmp@long_name	= "TOA CRF"
			var_tmp		= (/(FSNT_tmp-FSNTC_tmp)-(FLNT_tmp-FLNTC_tmp)/)

			var_tmp@units	= "W/m2"
			delete([/FSNT_tmp,FSNTC_tmp,FLNT_tmp,FLNTC_tmp/])
		end if
		if(var.eq."CrfSfc")then
			FSNS_tmp	= f_tmp[:]->FSNS
			FSNSC_tmp	= f_tmp[:]->FSNSC
			FLNS_tmp	= f_tmp[:]->FLNS
			FLNSC_tmp	= f_tmp[:]->FLNSC

			var_tmp		= FSNS_tmp
			var_tmp@long_name	= "Surface CRF"
			var_tmp		= (/(FSNS_tmp-FSNSC_tmp)-(FLNS_tmp-FLNSC_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FSNS_tmp,FSNSC_tmp,FLNS_tmp,FLNSC_tmp/])
		end if
		if(var.eq."ACRF")then
			FSNT_tmp	= f_tmp[:]->FSNT
			FSNTC_tmp	= f_tmp[:]->FSNTC
			FLNT_tmp	= f_tmp[:]->FLNT
			FLNTC_tmp	= f_tmp[:]->FLNTC

			FSNS_tmp	= f_tmp[:]->FSNS
			FSNSC_tmp	= f_tmp[:]->FSNSC
			FLNS_tmp	= f_tmp[:]->FLNS
			FLNSC_tmp	= f_tmp[:]->FLNSC

			var_tmp		= FSNT_tmp
			var_tmp@long_name	= "Atm CRF"
			var_tmp		= (/((FSNT_tmp-FSNTC_tmp)-(FLNT_tmp-FLNTC_tmp))-((FSNS_tmp-FSNSC_tmp)-(FLNS_tmp-FLNSC_tmp))/)
			var_tmp@units	= "W/m2"
			delete([/FSNT_tmp,FSNTC_tmp,FLNT_tmp,FLNTC_tmp,FSNS_tmp,FSNSC_tmp,FLNS_tmp,FLNSC_tmp/])
		end if

		if(var.eq."RadCooling_allsky")then
			FSNT_tmp	= f_tmp[:]->FSNT
			FLNT_tmp	= f_tmp[:]->FLNT
			FSNS_tmp	= f_tmp[:]->FSNS
			FLNS_tmp	= f_tmp[:]->FLNS

			var_tmp		= FSNT_tmp
			var_tmp@long_name	= "Radiative cooling from allsky"
			var_tmp		= (/(FSNT_tmp-FLNT_tmp)-(FSNS_tmp-FLNS_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FSNT_tmp,FLNT_tmp,FSNS_tmp,FLNS_tmp/])
		end if
		if(var.eq."RadCooling_clrsky")then
			FSNTC_tmp	= f_tmp[:]->FSNTC
			FLNTC_tmp	= f_tmp[:]->FLNTC
			FSNSC_tmp	= f_tmp[:]->FSNSC
			FLNSC_tmp	= f_tmp[:]->FLNSC

			var_tmp		= FSNTC_tmp
			var_tmp@long_name	= "Radiative cooling from clearsky"
			var_tmp		= (/(FSNTC_tmp-FLNTC_tmp)-(FSNSC_tmp-FLNSC_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FSNTC_tmp,FLNTC_tmp,FSNSC_tmp,FLNSC_tmp/])
		end if
		if(var.eq."SWCF_Toa")then
			FSNT_tmp	= f_tmp[:]->FSNT
			FSNTC_tmp	= f_tmp[:]->FSNTC

			var_tmp		= FSNT_tmp
			var_tmp@long_name	= "TOA SWCF"
			var_tmp		= (/FSNT_tmp-FSNTC_tmp/)
			var_tmp@units	= "W/m2"
			delete([/FSNT_tmp,FSNTC_tmp/])
		end if
		if(var.eq."SWCF_Sfc")then
			FSNS_tmp	= f_tmp[:]->FSNS
			FSNSC_tmp	= f_tmp[:]->FSNSC

			var_tmp		= FSNS_tmp
			var_tmp@long_name	= "Sfc SWCF"
			var_tmp		= (/FSNS_tmp-FSNSC_tmp/)
			var_tmp@units	= "W/m2"
			delete([/FSNS_tmp,FSNSC_tmp/])
		end if

		if(var.eq."LWCF_Toa")then
			FLNT_tmp	= f_tmp[:]->FLNT
			FLNTC_tmp	= f_tmp[:]->FLNTC

			var_tmp		= FLNT_tmp
			var_tmp@long_name	= "TOA LWCF"
			var_tmp		= (/(-1.)*(FLNT_tmp-FLNTC_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FLNT_tmp,FLNTC_tmp/])
		end if
		if(var.eq."LWCF_Sfc")then
			FLNS_tmp	= f_tmp[:]->FLNS
			FLNSC_tmp	= f_tmp[:]->FLNSC

			var_tmp		= FLNS_tmp
			var_tmp@long_name	= "Sfc LWCF"
			var_tmp		= (/(-1.)*(FLNS_tmp-FLNSC_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FLNS_tmp,FLNSC_tmp/])
		end if
		if(var.eq."NetSfcRad")then
			FSNS_tmp 	= f_tmp[:]->FSNS
			FLNS_tmp 	= f_tmp[:]->FLNS
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = FSNS_tmp-FLNS_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp/])
		end if
		if(var.eq."V925")then
			var_tmp = f_tmp[:]->V(:,{925},:,:)
		end if
	end if

	if(var.eq."CLDTOT".or.var.eq."CLDLOW".or.var.eq."CLDMED".or.var.eq."CLDHGH")then
		var_tmp 	= var_tmp*100.
	end if

	var_tmp2 = var_tmp

	if(isfilevar(f_tmp_2,"OCNFRAC"))then
		ocnfrac = f_tmp[:]->OCNFRAC(0,:,:)
		var_tmp2 = (/mask(var_tmp,ocnfrac.ge.0.1,True)/)
		var_tmp2 = (/mask(var_tmp,ocnfrac.ge.0.1,True)/)
	end if
	data_out = var_tmp2
	lat_out = var_tmp2&lat 
	lon_out = var_tmp2&lon
	time_out = var_tmp2&time
	delete([/var_tmp,var_tmp2/])

	
end if

return ([/data_out,lat_out,lon_out,time_out/])

end 
;========================================================================
;========================================================================

function get_var_3D(fdir:string,fapend:string,case:string,var:string,type:string)
begin

cp			= 1.00464e3 ; J/kg/K
Lv			= 2.501e6 ; J/kg
g			= 9.81

f_tmp			= addfile (fdir+case+fapend, "r")                         ; add file

if(type.eq."model")then
	if(isfilevar(f_tmp,var))then
		var_tmp = f_tmp->$var$(0,:,:,:)
	else
		if(var.eq."MSE")then
			T_tmp		= f_tmp->T(0,:,:,:)
			Q_tmp		= f_tmp->Q(0,:,:,:)
			Z3_tmp		= f_tmp->Z3(0,:,:,:)
						
			var_tmp = T_tmp
			var_tmp@long_name = "moist static energy (J/kg)"
			var_tmp = (/cp*T_tmp+g*Z3_tmp+Lv*Q_tmp/)
			var_tmp@units = "J/kg"
			delete([/T_tmp,Q_tmp,Z3_tmp/])
		end if
	end if
	var_tmp2 = var_tmp

;	if(isfilevar(f_tmp,"OCNFRAC"))then
;		ocnfrac = f_tmp->OCNFRAC(0,:,:)
;		var_tmp2 = (/mask(var_tmp,ocnfrac.ge.0.5,True)/)
;	end if

	data_out = var_tmp2
	lat_out = var_tmp2&lat 
	lon_out = var_tmp2&lon
	lev_out	= var_tmp2&lev
	delete([/var_tmp,var_tmp2/])

	
end if

if(type.eq."obs")then
	if(isfilevar(f_tmp,var))then
			var_tmp = f_tmp->$var$
			if(dimsizes(dimsizes(var_tmp)).ne.2)then
				delete(var_tmp)
				var_tmp = f_tmp->$var$(0,:,:)
			end if
	else
		if(var.eq."MSE")then
			T_tmp		= f_tmp->T(0,:,:,:)
			Q_tmp		= f_tmp->SHUM(0,:,:,:) ; g/kg
			Q_tmp		= Q_tmp/1000.
			Z3_tmp		= f_tmp->Z3(0,:,:,:)
						
			var_tmp = T_tmp
			var_tmp@long_name = "moist static energy (J/kg)"
			var_tmp = (/cp*T_tmp+g*Z3_tmp+Lv*Q_tmp/)
			var_tmp@units = "J/kg"
			delete([/T_tmp,Q_tmp,Z3_tmp/])
		end if
	end if

	data_out = var_tmp
	lat_out = var_tmp&lat
	lon_out = var_tmp&lon
	lev_out	= var_tmp&lev
	delete([/var_tmp/])
end if
	
return ([/data_out,lat_out,lon_out,lev_out/])

end 

;========================================================================
;========================================================================

function get_var_2D(fdir:string,fapend:string,case:string,var:string,type:string)
begin

 f_tmp     = addfile (fdir+case+fapend, "r")                         ; add file

if(type.eq."model")then
	if(isfilevar(f_tmp,var))then
		var_tmp = f_tmp->$var$(0,:,:)
		if(var.eq."TAUX".or.var.eq."TAUY")then
			var_tmp = var_tmp*(-1.)
		end if
		if(var.eq."PSL")then
			var_tmp	= var_tmp/100. ;change from Pa to hPa
		end if
	else
		if(var.eq."SST")then
			var_tmp1   = f_tmp->TEMP(0,0,:,:) ; time,lev,lat,lon
			var_tmp = PopLatLon(var_tmp1,"gx1v6","1x1d","bilin","da","100716")
			delete([/var_tmp1/])
		end if
		if(var.eq."PRECT")then
			PRECC_tmp   = f_tmp->PRECC(0,:,:)
			PRECL_tmp   = f_tmp->PRECL(0,:,:)
			
			var_tmp = PRECC_tmp
			var_tmp@long_name = "total precipitation rate"
			var_tmp   = PRECC_tmp+PRECL_tmp
			
			var_tmp = var_tmp*8.64e7
			var_tmp@units = "mm/day"
			
			delete([/PRECC_tmp,PRECL_tmp/])
		end if
		if(var.eq."LHFLX_revised")then
			PRECL_tmp	= f_tmp->PRECL(0,:,:)
			PRECC_tmp	= f_tmp->PRECC(0,:,:)
			PRECSC_tmp	= f_tmp->PRECSC(0,:,:)
			PRECSL_tmp 	= f_tmp->PRECSL(0,:,:)
			QFLX_tmp	= f_tmp->QFLX(0,:,:)
			
			Lv = 2.501e6                   ; J/kg
			Lf = 3.337e5                   ; J/kg
			LHFLX_new = (/(Lv+Lf)*QFLX_tmp - Lf*1.e3*(PRECC_tmp+PRECL_tmp-PRECSC_tmp-PRECSL_tmp)/)
			
			var_tmp = QFLX_tmp
			var_tmp@long_name = "revised LHFLX"
			var_tmp = (/LHFLX_new/)
			var_tmp@units = "W/m2"
			delete([/PRECL_tmp,PRECC_tmp,PRECSC_tmp,PRECSL_tmp,QFLX_tmp/])
		end if

		if(var.eq."NetSfcHeat")then
			FSNS_tmp 	= f_tmp->FSNS(0,:,:)
			FLNS_tmp 	= f_tmp->FLNS(0,:,:)
			LHFLX_tmp 	= f_tmp->LHFLX(0,:,:)
			SHFLX_tmp 	= f_tmp->SHFLX(0,:,:)

			PRECL_tmp	= f_tmp->PRECL(0,:,:)
			PRECC_tmp	= f_tmp->PRECC(0,:,:)
			PRECSC_tmp	= f_tmp->PRECSC(0,:,:)
			PRECSL_tmp 	= f_tmp->PRECSL(0,:,:)
			QFLX_tmp	= f_tmp->QFLX(0,:,:)
			
			Lv = 2.501e6                   ; J/kg
			Lf = 3.337e5                   ; J/kg
			LHFLX_new = (/(Lv+Lf)*QFLX_tmp - Lf*1.e3*(PRECC_tmp+PRECL_tmp-PRECSC_tmp-PRECSL_tmp)/)
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = FSNS_tmp-FLNS_tmp-LHFLX_new-SHFLX_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp,LHFLX_tmp,SHFLX_tmp/])
		end if
		if(var.eq."NetAtmHeat")then
			FSNS_tmp 	= f_tmp->FSNS(0,:,:)
			FLNS_tmp 	= f_tmp->FLNS(0,:,:)
			LHFLX_tmp 	= f_tmp->LHFLX(0,:,:)
			SHFLX_tmp 	= f_tmp->SHFLX(0,:,:)
			FSNT_tmp 	= f_tmp->FSNT(0,:,:)
			FLUT_tmp 	= f_tmp->FLNT(0,:,:)

			PRECL_tmp	= f_tmp->PRECL(0,:,:)
			PRECC_tmp	= f_tmp->PRECC(0,:,:)
			PRECSC_tmp	= f_tmp->PRECSC(0,:,:)
			PRECSL_tmp 	= f_tmp->PRECSL(0,:,:)
			QFLX_tmp	= f_tmp->QFLX(0,:,:)
			
			Lv = 2.501e6                   ; J/kg
			Lf = 3.337e5                   ; J/kg
			LHFLX_new = (/(Lv+Lf)*QFLX_tmp - Lf*1.e3*(PRECC_tmp+PRECL_tmp-PRECSC_tmp-PRECSL_tmp)/)
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net atm heat flux"
			var_tmp = (FSNT_tmp-FLUT_tmp)-(FSNS_tmp-FLNS_tmp-LHFLX_new-SHFLX_tmp)
			var_tmp@units = "W/m2"
			delete([/FSNT_tmp,FLUT_tmp,FSNS_tmp,FLNS_tmp,LHFLX_tmp,SHFLX_tmp/])
		end if
		if(var.eq."NetToaHeat")then
			FSNT_tmp 	= f_tmp->FSNT(0,:,:)
			FLUT_tmp 	= f_tmp->FLNT(0,:,:)

			var_tmp = FSNT_tmp
			var_tmp@long_name = "Net TOA heat flux"
			var_tmp = (/FSNT_tmp-FLUT_tmp/)
			var_tmp@units = "W/m2"
			delete([/FSNT_tmp,FLUT_tmp/])
		end if
		if(var.eq."CrfToa")then
			FSNT_tmp	= f_tmp->FSNT(0,:,:)
			FSNTC_tmp	= f_tmp->FSNTC(0,:,:)
			FLNT_tmp	= f_tmp->FLNT(0,:,:)
			FLNTC_tmp	= f_tmp->FLNTC(0,:,:)

			var_tmp		= FSNT_tmp
			var_tmp@long_name	= "TOA CRF"
			var_tmp		= (/(FSNT_tmp-FSNTC_tmp)-(FLNT_tmp-FLNTC_tmp)/)
;			var_tmp		= (/(FSNT_tmp-FSNTC_tmp)/)

			var_tmp@units	= "W/m2"
			delete([/FSNT_tmp,FSNTC_tmp,FLNT_tmp,FLNTC_tmp/])
		end if
		if(var.eq."CrfSfc")then
			FSNS_tmp	= f_tmp->FSNS(0,:,:)
			FSNSC_tmp	= f_tmp->FSNSC(0,:,:)
			FLNS_tmp	= f_tmp->FLNS(0,:,:)
			FLNSC_tmp	= f_tmp->FLNSC(0,:,:)

			var_tmp		= FSNS_tmp
			var_tmp@long_name	= "Surface CRF"
			var_tmp		= (/(FSNS_tmp-FSNSC_tmp)-(FLNS_tmp-FLNSC_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FSNS_tmp,FSNSC_tmp,FLNS_tmp,FLNSC_tmp/])
		end if
		if(var.eq."ACRF")then
			FSNT_tmp	= f_tmp->FSNT(0,:,:)
			FSNTC_tmp	= f_tmp->FSNTC(0,:,:)
			FLNT_tmp	= f_tmp->FLNT(0,:,:)
			FLNTC_tmp	= f_tmp->FLNTC(0,:,:)

			FSNS_tmp	= f_tmp->FSNS(0,:,:)
			FSNSC_tmp	= f_tmp->FSNSC(0,:,:)
			FLNS_tmp	= f_tmp->FLNS(0,:,:)
			FLNSC_tmp	= f_tmp->FLNSC(0,:,:)

			var_tmp		= FSNT_tmp
			var_tmp@long_name	= "Atm CRF"
			var_tmp		= (/((FSNT_tmp-FSNTC_tmp)-(FLNT_tmp-FLNTC_tmp))-((FSNS_tmp-FSNSC_tmp)-(FLNS_tmp-FLNSC_tmp))/)
			var_tmp@units	= "W/m2"
			delete([/FSNT_tmp,FSNTC_tmp,FLNT_tmp,FLNTC_tmp,FSNS_tmp,FSNSC_tmp,FLNS_tmp,FLNSC_tmp/])
		end if

		if(var.eq."RadCooling_allsky")then
			FSNT_tmp	= f_tmp->FSNT(0,:,:)
			FLNT_tmp	= f_tmp->FLNT(0,:,:)
			FSNS_tmp	= f_tmp->FSNS(0,:,:)
			FLNS_tmp	= f_tmp->FLNS(0,:,:)

			var_tmp		= FSNT_tmp
			var_tmp@long_name	= "Radiative cooling from allsky"
			var_tmp		= (/(FSNT_tmp-FLNT_tmp)-(FSNS_tmp-FLNS_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FSNT_tmp,FLNT_tmp,FSNS_tmp,FLNS_tmp/])
		end if
		if(var.eq."RadCooling_clrsky")then
			FSNTC_tmp	= f_tmp->FSNTC(0,:,:)
			FLNTC_tmp	= f_tmp->FLNTC(0,:,:)
			FSNSC_tmp	= f_tmp->FSNSC(0,:,:)
			FLNSC_tmp	= f_tmp->FLNSC(0,:,:)

			var_tmp		= FSNTC_tmp
			var_tmp@long_name	= "Radiative cooling from clearsky"
			var_tmp		= (/(FSNTC_tmp-FLNTC_tmp)-(FSNSC_tmp-FLNSC_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FSNTC_tmp,FLNTC_tmp,FSNSC_tmp,FLNSC_tmp/])
		end if
		if(var.eq."SWCF_Toa")then
			FSNT_tmp	= f_tmp->FSNT(0,:,:)
			FSNTC_tmp	= f_tmp->FSNTC(0,:,:)

			var_tmp		= FSNT_tmp
			var_tmp@long_name	= "TOA SWCF"
			var_tmp		= (/FSNT_tmp-FSNTC_tmp/)
			var_tmp@units	= "W/m2"
			delete([/FSNT_tmp,FSNTC_tmp/])
		end if
		if(var.eq."SWCF_Sfc")then
			FSNS_tmp	= f_tmp->FSNS(0,:,:)
			FSNSC_tmp	= f_tmp->FSNSC(0,:,:)

			var_tmp		= FSNS_tmp
			var_tmp@long_name	= "Sfc SWCF"
			var_tmp		= (/FSNS_tmp-FSNSC_tmp/)
			var_tmp@units	= "W/m2"
			delete([/FSNS_tmp,FSNSC_tmp/])
		end if

		if(var.eq."LWCF_Toa")then
			FLNT_tmp	= f_tmp->FLNT(0,:,:)
			FLNTC_tmp	= f_tmp->FLNTC(0,:,:)

			var_tmp		= FLNT_tmp
			var_tmp@long_name	= "TOA LWCF"
			var_tmp		= (/(-1.)*(FLNT_tmp-FLNTC_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FLNT_tmp,FLNTC_tmp/])
		end if
		if(var.eq."LWCF_Sfc")then
			FLNS_tmp	= f_tmp->FLNS(0,:,:)
			FLNSC_tmp	= f_tmp->FLNSC(0,:,:)

			var_tmp		= FLNS_tmp
			var_tmp@long_name	= "Sfc LWCF"
			var_tmp		= (/(-1.)*(FLNS_tmp-FLNSC_tmp)/)
			var_tmp@units	= "W/m2"
			delete([/FLNS_tmp,FLNSC_tmp/])
		end if
		if(var.eq."LH_SH")then
			SHFLX_tmp 	= f_tmp->SHFLX(0,:,:)

			PRECL_tmp	= f_tmp->PRECL(0,:,:)
			PRECC_tmp	= f_tmp->PRECC(0,:,:)
			PRECSC_tmp	= f_tmp->PRECSC(0,:,:)
			PRECSL_tmp 	= f_tmp->PRECSL(0,:,:)
			QFLX_tmp	= f_tmp->QFLX(0,:,:)
			
			Lv = 2.501e6                   ; J/kg
			Lf = 3.337e5                   ; J/kg
			LHFLX_new = (/(Lv+Lf)*QFLX_tmp - Lf*1.e3*(PRECC_tmp+PRECL_tmp-PRECSC_tmp-PRECSL_tmp)/)
			
			var_tmp = SHFLX_tmp
			var_tmp@long_name = "latent+sensible heat flux"
			var_tmp = (/LHFLX_new+SHFLX_tmp/)
			var_tmp@units = "W/m2"
			delete([/SHFLX_tmp,PRECL_tmp,PRECC_tmp,PRECSC_tmp,PRECSL_tmp,QFLX_tmp,LHFLX_new/])
		end if

		if(var.eq."NetSfcRad")then
			FSNS_tmp 	= f_tmp->FSNS(0,:,:)
			FLNS_tmp 	= f_tmp->FLNS(0,:,:)
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = FSNS_tmp-FLNS_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp/])
		end if
		if(var.eq."V925")then
			var_tmp = f_tmp->V(0,{925},:,:)
		end if
		if(var.eq."U925")then
			var_tmp = f_tmp->U(0,{925},:,:)
		end if

	end if

	if(var.eq."CLDTOT".or.var.eq."CLDLOW".or.var.eq."CLDMED".or.var.eq."CLDHGH")then
		var_tmp 	= var_tmp*100.
	end if

	var_tmp2 = var_tmp
	
;	if(var.ne."PS")then ; "surface pressure" cannot be allowed to be masked. so ignore masking it here.
;		if(isfilevar(f_tmp,"OCNFRAC"))then
;			ocnfrac = f_tmp->OCNFRAC(0,:,:)
;			var_tmp2 = (/mask(var_tmp,ocnfrac.ge.0.1,True)/)
;			var_tmp2 = (/mask(var_tmp,ocnfrac.ge.0.1,True)/)
;		end if
;	end if

	data_out = var_tmp2
	lat_out = var_tmp2&lat 
	lon_out = var_tmp2&lon
	delete([/var_tmp,var_tmp2/])

	
end if

if(type.eq."obs")then
	if(isfilevar(f_tmp,var))then
			var_tmp = f_tmp->$var$
			if(dimsizes(dimsizes(var_tmp)).ne.2)then
				delete(var_tmp)
				var_tmp = f_tmp->$var$(0,:,:)
			end if
	else
		if(var.eq."SST")then
			var_tmp = f_tmp->SST
		end if	
		if(var.eq."PRECT")then
			var_tmp = f_tmp->PRECT(0,:,:)
		end if
		if(var.eq."NetSfcHeat")then
			FSNS_tmp 	= f_tmp->FSNS(0,:,:)
			FLNS_tmp 	= f_tmp->FLNS(0,:,:)
			LHFLX_tmp 	= f_tmp->LHFLX(0,:,:)
			SHFLX_tmp 	= f_tmp->SHFLX(0,:,:)
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = (-1.0)*FSNS_tmp-FLNS_tmp-LHFLX_tmp-SHFLX_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp,LHFLX_tmp,SHFLX_tmp/])
		end if
		if(var.eq."NetAtmHeat")then
			FSNT_tmp 	= f_tmp->FSNT(0,:,:)
			FLUT_tmp 	= f_tmp->FLUT(0,:,:)
			FSNS_tmp 	= f_tmp->FSNS(0,:,:)
			FLNS_tmp 	= f_tmp->FLNS(0,:,:)
			LHFLX_tmp 	= f_tmp->LHFLX(0,:,:)
			SHFLX_tmp 	= f_tmp->SHFLX(0,:,:)
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = ((-1.0)*FSNT_tmp-FLUT_tmp)-((-1.)*FSNS_tmp-FLNS_tmp-LHFLX_tmp-SHFLX_tmp)
			var_tmp@units = "W/m2"
			delete([/FSNT_tmp,FLUT_tmp,FSNS_tmp,FLNS_tmp,LHFLX_tmp,SHFLX_tmp/])
		end if

		if(var.eq."NetToaHeat")then
			FSNT_tmp 	= f_tmp->FSNT(0,:,:)
			FLUT_tmp 	= f_tmp->FLUT(0,:,:)

			var_tmp = FSNT_tmp
			var_tmp@long_name = "Net TOA heat flux"
			var_tmp = (-1.)*FSNT_tmp-FLUT_tmp
			var_tmp@units = "W/m2"
			delete([/FSNT_tmp,FLUT_tmp/])
		end if

		if(var.eq."NetSfcRad")then
			FSNS_tmp 	= f_tmp->FSNS(0,:,:)
			FLNS_tmp 	= f_tmp->FLNS(0,:,:)
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = (-1.0)*FSNS_tmp-FLNS_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp/])
		end if
		if(var.eq."CrfToa")then
			FSNT_tmp	= f_tmp->FSUT(0,:,:)
			FSNTC_tmp	= f_tmp->FSUTC(0,:,:)
			FLNT_tmp	= f_tmp->FLUT(0,:,:)
			FLNTC_tmp	= f_tmp->FLUTC(0,:,:)

			var_tmp		= FSNT_tmp
			var_tmp@long_name	= "TOA CRF"
			var_tmp		= (/(-1.)*(FSNT_tmp-FSNTC_tmp)-(FLNT_tmp-FLNTC_tmp)/)
;			var_tmp		= (/(-1.)*(FSNT_tmp-FSNTC_tmp)/)

			var_tmp@units	= "W/m2"
			delete([/FSNT_tmp,FSNTC_tmp,FLNT_tmp,FLNTC_tmp/])
		end if

		if(var.eq."V925")then
			var_tmp = f_tmp->V(0,{925},:,:)
		end if
		if(var.eq."U925")then
			var_tmp = f_tmp->U(0,{925},:,:)
		end if

	end if

	if(var.eq."FSNS".or.var.eq."FSNT")then
		var_tmp = var_tmp*(-1.)
	end if

	data_out = var_tmp
	lat_out = var_tmp&lat
	lon_out = var_tmp&lon
	delete([/var_tmp/])
end if
	
return ([/data_out,lat_out,lon_out/])

end 
;========================================================================
;========================================================================

function get_var_time(fdir:string,fapend:string,case:string,var:string,type:string)
begin

 f_tmp     = addfile (fdir+case+fapend, "r")                         ; add file

if(type.eq."model")then

	if(isfilevar(f_tmp,var))then
		var_tmp = f_tmp->$var$
	else
		if(var.eq."SST")then
			var_tmp1   = f_tmp->TEMP(:,0,:,:) ; time,lev,lat,lon
			var_tmp = PopLatLon(var_tmp1,"gx1v6","1x1d","bilin","da","100716")
			delete([/var_tmp1/])
		end if
		if(var.eq."PRECT")then
			PRECC_tmp   = f_tmp->PRECC
			PRECL_tmp   = f_tmp->PRECL
			
			var_tmp = PRECC_tmp
			var_tmp@long_name = "total precipitation rate"
			var_tmp   = PRECC_tmp+PRECL_tmp
			
			var_tmp = var_tmp*8.64e7
			var_tmp@units = "mm/day"
			
			delete([/PRECC_tmp,PRECL_tmp/])
		end if
		if(var.eq."NetSfcHeat")then
			FSNS_tmp 	= f_tmp->FSNS
			FLNS_tmp 	= f_tmp->FLNS
			LHFLX_tmp 	= f_tmp->LHFLX
			SHFLX_tmp 	= f_tmp->SHFLX
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = FSNS_tmp-FLNS_tmp-LHFLX_tmp-SHFLX_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp,LHFLX_tmp,SHFLX_tmp/])
		end if
		if(var.eq."NetSfcRad")then
			FSNS_tmp 	= f_tmp->FSNS
			FLNS_tmp 	= f_tmp->FLNS
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = FSNS_tmp-FLNS_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp/])
		end if
		if(var.eq."V925")then
			var_tmp = f_tmp->V(:,{925},:,:)
		end if
	end if

	if(var.eq."CLDTOT".or.var.eq."CLDLOW".or.var.eq."CLDMED".or.var.eq."CLDHGH")then
		var_tmp 	= var_tmp*100.
	end if


data_out 	= var_tmp
lat_out 	= var_tmp&lat 
lon_out 	= var_tmp&lon
time_out	= var_tmp&time
delete([/var_tmp/])

	
end if

if(type.eq."obs")then

	if(isfilevar(f_tmp,var))then
		var_tmp = f_tmp->$var$
	else
		if(var.eq."SST")then
			var_tmp = f_tmp->SST
			var_tmp!0 = "time"
			var_tmp&time = fspan(0,11,12)
		end if	
		if(var.eq."PRECT")then
			var_tmp = f_tmp->PRECT
		end if
		if(var.eq."NetSfcHeat")then
			FSNS_tmp 	= f_tmp->FSNS
			FLNS_tmp 	= f_tmp->FLNS
			LHFLX_tmp 	= f_tmp->LHFLX
			SHFLX_tmp 	= f_tmp->SHFLX
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = (-1.)*FSNS_tmp-FLNS_tmp-LHFLX_tmp-SHFLX_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp,LHFLX_tmp,SHFLX_tmp/])
		end if
		if(var.eq."NetSfcRad")then
			FSNS_tmp 	= f_tmp->FSNS
			FLNS_tmp 	= f_tmp->FLNS
			
			var_tmp = FSNS_tmp
			var_tmp@long_name = "Net surface heat flux"
			var_tmp = (-1.)*FSNS_tmp-FLNS_tmp
			var_tmp@units = "W/m2"
			delete([/FSNS_tmp,FLNS_tmp/])
		end if
		if(var.eq."V925")then
			var_tmp = f_tmp->V(:,{925},:,:)
		end if
	end if
		
	if(var.eq."FSNS")then
		var_tmp = var_tmp*(-1.)
	end if
	
	data_out 	= var_tmp
	lat_out 	= var_tmp&lat
	lon_out 	= var_tmp&lon
	time_out 	= var_tmp&time
	delete([/var_tmp/])

end if
	
return ([/data_out,lat_out,lon_out,time_out/])

end 

;========================================================================
;========================================================================

function get_wind_2D(fdir:string,fapend:string,case:string,var:string,type:string)
begin

 f_tmp     = addfile (fdir+case+fapend, "r")                         ; add file

if(type.eq."model")then
	U_tmp = f_tmp->U(0,{925},:,:)
	V_tmp = f_tmp->V(0,{925},:,:)

	data_out1 = U_tmp
	data_out2 = V_tmp
	lat_out = U_tmp&lat 
	lon_out = U_tmp&lon

	delete([/U_tmp,V_tmp/])
end if

if(type.eq."obs")then
	U_tmp = f_tmp->U(0,{925},:,:)
	V_tmp = f_tmp->V(0,{925},:,:)

	data_out1 = U_tmp
	data_out2 = V_tmp
	lat_out = U_tmp&lat 
	lon_out = U_tmp&lon

	delete([/U_tmp,V_tmp/])
end if
	
return ([/data_out1,data_out2,lat_out,lon_out/])

end 
;========================================================================
;========================================================================

function get_stress_2D(fdir:string,fapend:string,case:string,var:string,type:string)
begin

 f_tmp     = addfile (fdir+case+fapend, "r")                         ; add file

if(type.eq."model")then
	taux = f_tmp->TAUX(0,:,:)
	tauy = f_tmp->TAUY(0,:,:)

	tx	= taux 
	ty	= tauy
	if(isfilevar(f_tmp,"OCNFRAC"))then
		ocnfrac = f_tmp->OCNFRAC(0,:,:)
		tx = (/-1.0*mask(taux,ocnfrac.ge.0.1,True)/)
		ty = (/-1.0*mask(tauy,ocnfrac.ge.0.1,True)/)
	end if
	taux	= tx
	tauy	= ty
	delete(tx)
	delete(ty)

	data_out1 = taux
	data_out2 = tauy
	lat_out = taux&lat 
	lon_out = tauy&lon

	delete([/taux,tauy/])
end if

if(type.eq."obs")then
	U_tmp = f_tmp->TAUX(0,:,:)
	V_tmp = f_tmp->TAUY(0,:,:)

	data_out1 = U_tmp
	data_out2 = V_tmp
	lat_out = U_tmp&lat 
	lon_out = U_tmp&lon

	delete([/U_tmp,V_tmp/])
end if
	
return ([/data_out1,data_out2,lat_out,lon_out/])

end 


;========================================================================
;========================================================================

; 2017-10-26
function get_var_ocn_3D_depth(fdir:string,fapend:string,case:string,var:string,type:string)
begin

 setfileoption("nc","MissingToFillValue",False) 
 f_tmp     = addfile (fdir+case+fapend, "r")                         ; add file

if(type.eq."model")then
	if(isfilevar(f_tmp,var))then
		var_tmp1 	=	f_tmp->$var$(0,:,:,:)
		var_tmp		=	PopLatLon(var_tmp1,"gx1v6","1x1d","bilin","da","100716") 
	end if
	
	if(var.eq."UVEL".or.var.eq."VVEL".or.var.eq."WVEL")then
		var_tmp 	= 	var_tmp*1e-2	; change unit from "cm/s" to "m/s"
	end if

	; get depth variable
	z_t = f_tmp->z_t
	z_t = z_t*1e-2
	z_t@longname = "Depth"
	z_t@units = "m"

	data_out = var_tmp
	lat_out = var_tmp&lat 
	lon_out = var_tmp&lon
	depth_out = z_t

	delete([/var_tmp,z_t/])
	
end if

; 2018-01-06
if(type.eq."obs")then
	if(isfilevar(f_tmp,var))then
			var_tmp_db = f_tmp->$var$(0,:,:,:)
	end if
	var_tmp			= tofloat(var_tmp_db)
	copy_VarCoords(var_tmp_db,var_tmp)

	; change coordinates name mannually
	var_tmp!0		= "depth"
	var_tmp!1		= "lat"
	var_tmp!2		= "lon"
	var_tmp@_FillValue	= -1e20

	
	data_out 		= var_tmp
	lat_out 		= var_tmp&lat
	lon_out 		= var_tmp&lon
	depth_out 		= var_tmp&depth
	delete([/var_tmp/])
end if
	
return ([/data_out,lat_out,lon_out,depth_out/])

end 
;========================================================================
;========================================================================

;2017-10-27
function get_var_ocn_2D(fdir:string,fapend:string,case:string,var:string,type:string)
begin

 f_tmp     = addfile (fdir+case+fapend, "r")                         ; add file

if(type.eq."model")then
	if(isfilevar(f_tmp,var))then
		var_tmp1 	=	f_tmp->$var$(0,:,:)
		var_tmp		=	PopLatLon(var_tmp1,"gx1v6","1x1d","bilin","da","100716") 
	end if

	if(var.eq."HMXL")then
		var_tmp = var_tmp/100.		;change unit from "cm" to "m".
		var_tmp@units = "m"
	end if

	data_out = var_tmp
	lat_out = var_tmp&lat 
	lon_out = var_tmp&lon

	delete([/var_tmp/])
	
end if

;if(type.eq."obs")then
;	if(isfilevar(f_tmp,var))then
;			var_tmp = f_tmp->$var$
;			if(dimsizes(dimsizes(var_tmp)).ne.2)then
;				delete(var_tmp)
;				var_tmp = f_tmp->$var$(0,:,:)
;			end if
;	else
;		if(var.eq."SST")then
;			var_tmp = f_tmp->SST
;		end if	
;	end if
;
;	data_out = var_tmp
;	lat_out = var_tmp&lat
;	lon_out = var_tmp&lon
;	delete([/var_tmp/])
;end if
	
return ([/data_out,lat_out,lon_out/])

end 
;========================================================================
;========================================================================

function get_var_ocn_TimeLatLon_addfiles(Dir:string,fdir_append:string,case:string,var:string,type:string)
begin

fils 			= systemfunc ("ls "+Dir+case+fdir_append+"*pop.h.0*.nc") ; file paths
;fils 			= systemfunc ("ls "+Dir+case+fdir_append+"*pop.h.0001*.nc") ; file paths
;fils 			= systemfunc ("ls "+Dir+case+fdir_append+"*pop.h.000[123456]*.nc") ; file paths
;fils 			= systemfunc ("ls "+Dir+case+fdir_append+"*pop.h.000[12]*.nc") ; file paths

f_tmp     		= addfiles (fils, "r")                         ; add file

fils_2			= systemfunc ("ls "+Dir+case+fdir_append+"*0001-01.nc") ; file paths
f_tmp_2			= addfile(fils_2,"r")

ListSetType (f_tmp, "cat")

if(type.eq."model")then
	if(isfilevar(f_tmp_2,var))then
		var_tmp1 	=	f_tmp[:]->$var$
		var_tmp		=	PopLatLon(var_tmp1,"gx1v6","1x1d","bilin","da","100716") 
;		var_tmp		=	PopLatLon(var_tmp1,"gx1v6","fv1.9x2.5","bilin","da","090206")
	end if

	if(var.eq."HMXL")then
		var_tmp = var_tmp/100.		;change unit from "cm" to "m".
		var_tmp@units = "m"
	end if

	data_out = var_tmp
	lat_out = var_tmp&lat 
	lon_out = var_tmp&lon
	time_out = var_tmp&time
	delete([/var_tmp/])

end if

return ([/data_out,lat_out,lon_out,time_out/])

end 

; 2017-12-17
function get_var_ocn_TimeDepthLatLon_addfiles(Dir:string,fdir_append:string,case:string,var:string,type:string)
begin

fils 			= systemfunc ("ls "+Dir+case+fdir_append+"*pop.h.0*.nc") ; file paths
;fils 			= systemfunc ("ls "+Dir+case+fdir_append+"*pop.h.0001*.nc") ; file paths
;fils 			= systemfunc ("ls "+Dir+case+fdir_append+"*pop.h.000[123456]*.nc") ; file paths
;fils 			= systemfunc ("ls "+Dir+case+fdir_append+"*pop.h.000[12]*.nc") ; file paths

f_tmp     		= addfiles (fils, "r")                         ; add file

fils_2			= systemfunc ("ls "+Dir+case+fdir_append+"*0001-01.nc") ; file paths
f_tmp_2			= addfile(fils_2,"r")

ListSetType (f_tmp, "cat")

if(type.eq."model")then
	if(isfilevar(f_tmp_2,var))then
		var_tmp1 	=	f_tmp[:]->$var$
		var_tmp		=	PopLatLon(var_tmp1,"gx1v6","1x1d","bilin","da","100716") 
;		var_tmp		=	PopLatLon(var_tmp1,"gx1v6","fv1.9x2.5","bilin","da","090206") 
	end if

	if(var.eq."UVEL".or.var.eq."VVEL".or.var.eq."WVEL")then
		var_tmp 	= 	var_tmp*1e-2	; change unit from "cm/s" to "m/s"
	end if

	; get depth variable
	z_t = f_tmp_2->z_t
	z_t = z_t*1e-2
	z_t@longname = "Depth"
	z_t@units = "m"

	data_out = var_tmp
	lat_out = var_tmp&lat 
	lon_out = var_tmp&lon
	time_out = var_tmp&time
	depth_out = z_t
	delete([/var_tmp/])

end if
delete([/fils,f_tmp/])

return ([/data_out,lat_out,lon_out,time_out,depth_out/])

end 


