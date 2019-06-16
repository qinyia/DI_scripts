
function lat_lon_interpolate_2vars (a[*][*]:numeric,lat1[*]:numeric,lon1[*]:numeric,b[*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric)

begin

nlat1 = dimsizes(lat1)
nlat2 = dimsizes(lat2)

nlon1 = dimsizes(lon1)
nlon2 = dimsizes(lon2)


; =======find the minlatimum lat number
minlat = nlat1
minlon = nlon1
pointer = 1

if(minlat.ge.nlat2)then
	minlat = nlat2
	minlon = nlon2
	pointer = 2
end if

if(minlat.eq.nlat1)then
	lat_interp = lat1
	lon_interp = lon1
	copy_VarCoords(lat1,lat_interp)
	copy_VarCoords(lon1,lon_interp)
end if
if(minlat.eq.nlat2)then
	lat_interp = lat2
	lon_interp = lon2
	copy_VarCoords(lat2,lat_interp)
	copy_VarCoords(lon2,lon_interp)
end if


;==============interpolate data
plat = minlat
plon = minlon


ydata = new((/2,plat,plon/),"float")
ydata(0,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   
ydata(1,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)

  ydata!0 = "case"
  ydata&case = (/0,1/)
  ydata!1 = "lat"
  ydata&lat = lat_interp
  ydata!2 = "lon"
  ydata&lon = lon_interp

return(ydata)

end


function lat_interpolate (a[*]:numeric,lat1[*]:numeric,b[*]:numeric,lat2[*]:numeric,c[*]:numeric,lat3[*]:numeric)

begin

nlat1 = dimsizes(lat1)
nlat2 = dimsizes(lat2)
nlat3 = dimsizes(lat3)

; =======find the minlatimum lat number
minlat = nlat1
if(minlat.ge.nlat2)then
	minlat = nlat2
end if
if(minlat.ge.nlat3)then
	minlat = nlat3
end if

if(minlat.eq.nlat1)then
	lat_interp = lat1
	copy_VarCoords(lat1,lat_interp)
end if
if(minlat.eq.nlat2)then
	lat_interp = lat2
	copy_VarCoords(lat2,lat_interp)
end if
if(minlat.eq.nlat3)then
	lat_interp = lat3
	copy_VarCoords(lat3,lat_interp)
end if


;==============interpolate data
plat = minlat


ydata = new((/3,plat/),"float")
ydata(0,:) = linint1 (lat1,a,False,lat_interp,0)   
ydata(1,:) = linint1 (lat2,b,False,lat_interp,0)
ydata(2,:) = linint1 (lat3,c,False,lat_interp,0)

  ydata!0 = "line"
  ydata&line = (/0,1,2/)
  ydata!1 = "lat"
  ydata&lat = lat_interp

return(ydata)

end

; 2017-10-31
function lat_lon_interpolate_spcgrid (a[*][*]:numeric,lat1[*]:numeric,lon1[*]:numeric,\
			              b[*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric,\
				      c[*][*]:numeric,lat3[*]:numeric,lon3[*]:numeric,\
				      lat_spc[*]:numeric,lon_spc[*]:numeric)

begin

nlat  = dimsizes(lat_spc)
nlon  = dimsizes(lon_spc)

;==============interpolate data
plat = nlat
plon = nlon

lon_interp		= lon_spc
lat_interp		= lat_spc


ydata = new((/3,plat,plon/),"float")
ydata(0,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   
ydata(1,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)
ydata(2,:,:) = linint2 (lon3,lat3,c,False,lon_interp,lat_interp,0)

  ydata!0 = "case"
  ydata&case = (/0,1,2/)
  ydata!1 = "lat"
  ydata&lat = lat_interp
  ydata!2 = "lon"
  ydata&lon = lon_interp

return(ydata)

end

; 2017-11-01
function time_lat_lon_interpolate_2vars (ntime[*]:numeric,\
					a[*][*][*]:numeric,lat1[*]:numeric,lon1[*]:numeric,\
					b[*][*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric)

begin

nlat1 = dimsizes(lat1)
nlat2 = dimsizes(lat2)

nlon1 = dimsizes(lon1)
nlon2 = dimsizes(lon2)


; =======find the minlatimum lat number
minlat = nlat1
minlon = nlon1
pointer = 1

if(minlat.ge.nlat2)then
	minlat = nlat2
	minlon = nlon2
	pointer = 2
end if

if(minlat.eq.nlat1)then
	lat_interp = lat1
	lon_interp = lon1
	copy_VarCoords(lat1,lat_interp)
	copy_VarCoords(lon1,lon_interp)
end if
if(minlat.eq.nlat2)then
	lat_interp = lat2
	lon_interp = lon2
	copy_VarCoords(lat2,lat_interp)
	copy_VarCoords(lon2,lon_interp)
end if

;==============interpolate data
plat = minlat
plon = minlon


ydata = new((/2,ntime,plat,plon/),"float")


ydata(0,:,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   
ydata(1,:,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)

  ydata!0 = "case"
  ydata&case = (/0,1/)
  ydata!1 = "time"
  ydata&time = ispan(0,ntime-1,1)
  ydata!2 = "lat"
  ydata&lat = lat_interp
  ydata!3 = "lon"
  ydata&lon = lon_interp

return(ydata)

end

; 2017-12-17
function TimeDepthLatLon_interpolate_2vars (ntime[*]:numeric,ndepth[*]:numeric,\
					a[*][*][*][*]:numeric,lat1[*]:numeric,lon1[*]:numeric,\
					b[*][*][*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric)

begin

nlat1 = dimsizes(lat1)
nlat2 = dimsizes(lat2)

nlon1 = dimsizes(lon1)
nlon2 = dimsizes(lon2)


; =======find the minlatimum lat number
minlat = nlat1
minlon = nlon1
pointer = 1

if(minlat.ge.nlat2)then
	minlat = nlat2
	minlon = nlon2
	pointer = 2
end if

if(minlat.eq.nlat1)then
	lat_interp = lat1
	lon_interp = lon1
	copy_VarCoords(lat1,lat_interp)
	copy_VarCoords(lon1,lon_interp)
end if
if(minlat.eq.nlat2)then
	lat_interp = lat2
	lon_interp = lon2
	copy_VarCoords(lat2,lat_interp)
	copy_VarCoords(lon2,lon_interp)
end if

;==============interpolate data
plat = minlat
plon = minlon


ydata = new((/2,ntime,ndepth,plat,plon/),"float")


ydata(0,:,:,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   
print("ok 1")
ydata(1,:,:,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)
print("ok 2")

  ydata!0 = "case"
  ydata&case = (/0,1/)
  ydata!1 = "time"
  ydata&time = ispan(0,ntime-1,1)
  ydata!2 = "depth"
  ydata&depth = a&z_t
  ydata!3 = "lat"
  ydata&lat = lat_interp
  ydata!4 = "lon"
  ydata&lon = lon_interp

return(ydata)

end



function lat_lon_interpolate_3vars (a[*][*]:numeric,lat1[*]:numeric,lon1[*]:numeric,\
		                    b[*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric,\
				    c[*][*]:numeric,lat3[*]:numeric,lon3[*]:numeric)

begin

nlat1 = dimsizes(lat1)
nlat2 = dimsizes(lat2)
nlat3 = dimsizes(lat3)

nlon1 = dimsizes(lon1)
nlon2 = dimsizes(lon2)
nlon3 = dimsizes(lon3)


; =======find the minlatimum lat number
minlat = nlat1
minlon = nlon1
pointer = 1

if(minlat.ge.nlat2)then
	minlat = nlat2
	minlon = nlon2
	pointer = 2
end if
if(minlat.ge.nlat3)then
	minlat = nlat3
	minlon = nlon3
	pointer = 3
end if

if(minlat.eq.nlat1)then
	lat_interp = lat1
	lon_interp = lon1
	copy_VarCoords(lat1,lat_interp)
	copy_VarCoords(lon1,lon_interp)
end if
if(minlat.eq.nlat2)then
	lat_interp = lat2
	lon_interp = lon2
	copy_VarCoords(lat2,lat_interp)
	copy_VarCoords(lon2,lon_interp)
end if
if(minlat.eq.nlat3)then
	lat_interp = lat3
	lon_interp = lon3
	copy_VarCoords(lat3,lat_interp)
	copy_VarCoords(lon3,lon_interp)
end if



;==============interpolate data
plat = minlat
plon = minlon


ydata = new((/3,plat,plon/),"float")
ydata(0,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   

ydata(1,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)
ydata(2,:,:) = linint2 (lon3,lat3,c,False,lon_interp,lat_interp,0)

  ydata!0 = "case"
  ydata&case = (/0,1,2/)
  ydata!1 = "lat"
  ydata&lat = lat_interp
  ydata!2 = "lon"
  ydata&lon = lon_interp

return(ydata)

end

; 2017-12-12
function lat_lon_interpolate_3vars_specgrid (a[*][*]:numeric,lat1[*]:numeric,lon1[*]:numeric,\
		                   	     b[*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric,\
				   	     c[*][*]:numeric,lat3[*]:numeric,lon3[*]:numeric,\
					     lat_spc[*]:numeric,lon_spc[*]:numeric)

begin

;==============interpolate data
plat = dimsizes(lat_spc)
plon = dimsizes(lon_spc)

lat_interp = lat_spc
lon_interp = lon_spc

ydata = new((/3,plat,plon/),"float")
ydata(0,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   
ydata(1,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)
ydata(2,:,:) = linint2 (lon3,lat3,c,False,lon_interp,lat_interp,0)

  ydata!0 = "case"
  ydata&case = (/0,1,2/)
  ydata!1 = "lat"
  ydata&lat = lat_interp
  ydata!2 = "lon"
  ydata&lon = lon_interp

return(ydata)

end


;2017-11-02
function time_lat_lon_interpolate_3vars_specgrid (ntime[*]:numeric,\
						  a[*][*][*]:numeric,lat1[*]:numeric,lon1[*]:numeric,\
						  b[*][*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric,\
						  c[*][*][*]:numeric,lat3[*]:numeric,lon3[*]:numeric,\
						  lat_spc[*]:numeric,lon_spc[*]:numeric)

begin

;==============interpolate data
plat = dimsizes(lat_spc) 
plon = dimsizes(lon_spc)

lat_interp = lat_spc
lon_interp = lon_spc

ydata = new((/3,ntime,plat,plon/),"float")


ydata(0,:,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   
ydata(1,:,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)
ydata(2,:,:,:) = linint2 (lon3,lat3,c,False,lon_interp,lat_interp,0)

  ydata!0 = "case"
  ydata&case = (/0,1,2/)
  ydata!1 = "time"
  ydata&time = ispan(0,ntime-1,1)
  ydata!2 = "lat"
  ydata&lat = lat_interp
  ydata!3 = "lon"
  ydata&lon = lon_interp

return(ydata)

end



function time_lat_lon_interpolate_3vars (ntime[*]:numeric,\
					a[*][*][*]:numeric,lat1[*]:numeric,lon1[*]:numeric,\
					b[*][*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric,\
					c[*][*][*]:numeric,lat3[*]:numeric,lon3[*]:numeric)

begin

nlat1 = dimsizes(lat1)
nlat2 = dimsizes(lat2)
nlat3 = dimsizes(lat3)

nlon1 = dimsizes(lon1)
nlon2 = dimsizes(lon2)
nlon3 = dimsizes(lon3)


; =======find the minlatimum lat number
minlat = nlat1
minlon = nlon1
pointer = 1

if(minlat.ge.nlat2)then
	minlat = nlat2
	minlon = nlon2
	pointer = 2
end if
if(minlat.ge.nlat3)then
	minlat = nlat3
	minlon = nlon3
	pointer = 3
end if

if(minlat.eq.nlat1)then
	lat_interp = lat1
	lon_interp = lon1
	copy_VarCoords(lat1,lat_interp)
	copy_VarCoords(lon1,lon_interp)
end if
if(minlat.eq.nlat2)then
	lat_interp = lat2
	lon_interp = lon2
	copy_VarCoords(lat2,lat_interp)
	copy_VarCoords(lon2,lon_interp)
end if
if(minlat.eq.nlat3)then
	lat_interp = lat3
	lon_interp = lon3
	copy_VarCoords(lat3,lat_interp)
	copy_VarCoords(lon3,lon_interp)
end if



;==============interpolate data
plat = minlat
plon = minlon


ydata = new((/3,ntime,plat,plon/),"float")


ydata(0,:,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   
ydata(1,:,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)
ydata(2,:,:,:) = linint2 (lon3,lat3,c,False,lon_interp,lat_interp,0)

  ydata!0 = "case"
  ydata&case = (/0,1,2/)
  ydata!1 = "time"
  ydata&time = ispan(0,ntime-1,1)
  ydata!2 = "lat"
  ydata&lat = lat_interp
  ydata!3 = "lon"
  ydata&lon = lon_interp

return(ydata)

end

function lev_lat_lon_interpolate (a[*][*][*]:numeric,lev[*]:numeric,lat1[*]:numeric,lon1[*]:numeric,\
				  b[*][*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric,\
				  c[*][*][*]:numeric,lat3[*]:numeric,lon3[*]:numeric)

begin

nlat1 = dimsizes(lat1)
nlat2 = dimsizes(lat2)
nlat3 = dimsizes(lat3)

nlon1 = dimsizes(lon1)
nlon2 = dimsizes(lon2)
nlon3 = dimsizes(lon3)

nlev = dimsizes(lev)


; =======find the minlatimum lat number
minlat = nlat1
pointer = 1

if(minlat.ge.nlat2)then
	minlat = nlat2
	minlon = nlon2
	pointer = 2
end if
if(minlat.ge.nlat3)then
	minlat = nlat3
	minlon = nlon3
	pointer = 3
end if

if(minlat.eq.nlat1)then
	lat_interp = lat1
	lon_interp = lon1
	copy_VarCoords(lat1,lat_interp)
	copy_VarCoords(lon1,lon_interp)
end if
if(minlat.eq.nlat2)then
	lat_interp = lat2
	lon_interp = lon2
	copy_VarCoords(lat2,lat_interp)
	copy_VarCoords(lon2,lon_interp)
end if
if(minlat.eq.nlat3)then
	lat_interp = lat3
	lon_interp = lon3
	copy_VarCoords(lat3,lat_interp)
	copy_VarCoords(lon3,lon_interp)
end if



;==============interpolate data
plat = minlat
plon = minlon


ydata = new((/3,nlev,plat,plon/),"float")

ydata(0,:,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   
ydata(1,:,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)
ydata(2,:,:,:) = linint2 (lon3,lat3,c,False,lon_interp,lat_interp,0)

  ydata!0 = "case"
  ydata&case = (/0,1,2/)
  ydata!1 = "lev"
  ydata&lev = lev
  ydata!2 = "lat"
  ydata&lat = lat_interp
  ydata!3 = "lon"
  ydata&lon = lon_interp

return(ydata)

end

;2017-10-26
function lev_lat_lon_interpolate_2vars (a[*][*][*]:numeric,lev[*]:numeric,lat1[*]:numeric,lon1[*]:numeric,\
					b[*][*][*]:numeric,lat2[*]:numeric,lon2[*]:numeric)

begin

nlat1 = dimsizes(lat1)
nlat2 = dimsizes(lat2)

nlon1 = dimsizes(lon1)
nlon2 = dimsizes(lon2)

nlev = dimsizes(lev)


; =======find the minlatimum lat number
minlat = nlat1
minlon = nlon1
pointer = 1

if(minlat.ge.nlat2)then
	minlat = nlat2
	minlon = nlon2
	pointer = 2
end if

if(minlat.eq.nlat1)then
	lat_interp = lat1
	lon_interp = lon1
	copy_VarCoords(lat1,lat_interp)
	copy_VarCoords(lon1,lon_interp)
end if
if(minlat.eq.nlat2)then
	lat_interp = lat2
	lon_interp = lon2
	copy_VarCoords(lat2,lat_interp)
	copy_VarCoords(lon2,lon_interp)
end if


;==============interpolate data
plat = minlat
plon = minlon


ydata = new((/2,nlev,plat,plon/),"float")

ydata(0,:,:,:) = linint2 (lon1,lat1,a,False,lon_interp,lat_interp,0)   
ydata(1,:,:,:) = linint2 (lon2,lat2,b,False,lon_interp,lat_interp,0)

  ydata!0 = "case"
  ydata&case = (/0,1/)
  ydata!1 = "lev"
  ydata&lev = lev
  ydata!2 = "lat"
  ydata&lat = lat_interp
  ydata!3 = "lon"
  ydata&lon = lon_interp

return(ydata)

end
