
dir1				= "/home/lyl/WORK3/qinyi/double-ITCZ/data/obs/"
dir2				= "/home/lyl/WORK3/qinyi/double-ITCZ/data/atm/"

case				= (/"JRA25","BC5_f19g16","BC5_f19g16_mac2"/)
append				= "_ANN_climo_NetAtmHeat_regrid.nc"
dir				= (/dir1,dir2,dir2/)

nlat				= 180
nlon				= 360
NetAtmHeat			= new((/dimsizes(case),nlat,nlon/),double)
NetTOAHeat			= new((/dimsizes(case),nlat,nlon/),double)
NetSfcHeat			= new((/dimsizes(case),nlat,nlon/),double)

do icase=0,dimsizes(case)-1
ftmp				= addfile(dir(icase)+case(icase)+append,"r")

NetAtmHeat(icase,:,:)		= ftmp->NetAtmHeat
NetTOAHeat(icase,:,:)		= ftmp->NetTOAHeat
NetSfcHeat(icase,:,:)		= ftmp->NetSfcHeat

end do

printVarSummary(NetAtmHeat)

; weighting of latitude 

lat_here			= NetAtmHeat&lat
rad				= 4.0*atan(1.0)/180.
clat				= cos(lat_here*rad)
clat!0				= "lat"
clat&lat			= lat_here

latS1				= -90.
latE1				= 0.

latS2				= 0.
latE2				= 90.

lonS				= 0.
lonE				= 360.

; global mean
NetAtmHeat_GLB			= wgt_areaave_Wrap(NetAtmHeat(:,:,{lonS:lonE}),clat,1.0,1) 
NetTOAHeat_GLB			= wgt_areaave_Wrap(NetTOAHeat(:,:,{lonS:lonE}),clat,1.0,1) 
NetSfcHeat_GLB			= wgt_areaave_Wrap(NetSfcHeat(:,:,{lonS:lonE}),clat,1.0,1) 


; area mean over SH
NetAtmHeat_SH			= wgt_areaave_Wrap(NetAtmHeat(:,{latS1:latE1},{lonS:lonE}),clat({latS1:latE1}),1.0,1) 
NetTOAHeat_SH			= wgt_areaave_Wrap(NetTOAHeat(:,{latS1:latE1},{lonS:lonE}),clat({latS1:latE1}),1.0,1) 
NetSfcHeat_SH			= wgt_areaave_Wrap(NetSfcHeat(:,{latS1:latE1},{lonS:lonE}),clat({latS1:latE1}),1.0,1) 

; area mean over NH
NetAtmHeat_NH			= wgt_areaave_Wrap(NetAtmHeat(:,{latS2:latE2},{lonS:lonE}),clat({latS2:latE2}),1.0,1) 
NetTOAHeat_NH			= wgt_areaave_Wrap(NetTOAHeat(:,{latS2:latE2},{lonS:lonE}),clat({latS2:latE2}),1.0,1) 
NetSfcHeat_NH			= wgt_areaave_Wrap(NetSfcHeat(:,{latS2:latE2},{lonS:lonE}),clat({latS2:latE2}),1.0,1) 

; difference between SH and NH
NetAtmHeat_asy			= NetAtmHeat_SH-NetAtmHeat_NH
NetTOAHeat_asy			= NetTOAHeat_SH-NetTOAHeat_NH
NetSfcHeat_asy			= NetSfcHeat_SH-NetSfcHeat_NH
copy_VarCoords(NetAtmHeat_NH,NetAtmHeat_asy)
copy_VarCoords(NetTOAHeat_NH,NetTOAHeat_asy)
copy_VarCoords(NetSfcHeat_NH,NetSfcHeat_asy)

print("----------------------------------------JRA25--------------------------Park-RH------------------------Gauss-PDF---")
print("NetAtmHeat_GLB=			"+NetAtmHeat_GLB(0)+"		"+NetAtmHeat_GLB(1)+"		"+NetAtmHeat_GLB(2))
print("NetTOAHeat_GLB=			"+NetTOAHeat_GLB(0)+"		"+NetTOAHeat_GLB(1)+"		"+NetTOAHeat_GLB(2))
print("NetSfcHeat_GLB=			"+NetSfcHeat_GLB(0)+"		"+NetSfcHeat_GLB(1)+"		"+NetSfcHeat_GLB(2))

print("----------------------------------------JRA25--------------------------Park-RH------------------------Gauss-PDF---")
print("NetAtmHeat_SH=			"+NetAtmHeat_SH(0)+"		"+NetAtmHeat_SH(1)+"		"+NetAtmHeat_SH(2))
print("NetTOAHeat_SH=			"+NetTOAHeat_SH(0)+"		"+NetTOAHeat_SH(1)+"		"+NetTOAHeat_SH(2))
print("NetSfcHeat_SH=			"+NetSfcHeat_SH(0)+"		"+NetSfcHeat_SH(1)+"		"+NetSfcHeat_SH(2))

print("----------------------------------------JRA25--------------------------Park-RH------------------------Gauss-PDF---")
print("NetAtmHeat_NH=			"+NetAtmHeat_NH(0)+"		"+NetAtmHeat_NH(1)+"		"+NetAtmHeat_NH(2))
print("NetTOAHeat_NH=			"+NetTOAHeat_NH(0)+"		"+NetTOAHeat_NH(1)+"		"+NetTOAHeat_NH(2))
print("NetSfcHeat_NH=			"+NetSfcHeat_NH(0)+"		"+NetSfcHeat_NH(1)+"		"+NetSfcHeat_NH(2))

print("----------------------------------------JRA25--------------------------Park-RH------------------------Gauss-PDF---")
print("NetAtmHeat_asy=			"+NetAtmHeat_asy(0)+"		"+NetAtmHeat_asy(1)+"		"+NetAtmHeat_asy(2))
print("NetTOAHeat_asy=			"+NetTOAHeat_asy(0)+"		"+NetTOAHeat_asy(1)+"		"+NetTOAHeat_asy(2))
print("NetSfcHeat_asy=			"+NetSfcHeat_asy(0)+"		"+NetSfcHeat_asy(1)+"		"+NetSfcHeat_asy(2))


