
dir1				= "/home/lyl/WORK3/qinyi/double-ITCZ/data/obs/"
dir2				= "/home/lyl/WORK3/qinyi/double-ITCZ/data/atm/"

case				= (/"obs","BC5_f19g16","BC5_f19g16_mac2"/)
append				= "_glbvar_for_glbmean_regrid.nc"
dir				= (/dir1,dir2,dir2/)

nlat				= 180
nlon				= 360
vars				= (/"NetTOAHeat","FSNT","FLNT","SWCF","LWCF","CLDLOW","CLDMED","CLDHGH","CLDTOT"/)
; it seems that transition symbol cannot work now.
;vars_units			= (/"W/m~S~2~N~","W/m~S~2~N~","W/m~S~2~N~","W/m~S~2~N~","W/m~S~2~N~","%","%","%","%"/)
vars_units			= (/"W/m2","W/m2","W/m2","W/m2","W/m2","%","%","%","%"/)
print(vars_units)

VarsData			= new((/dimsizes(vars),dimsizes(case),nlat,nlon/),double)
VarsData!0			= "var"
VarsData&var			= ispan(0,dimsizes(vars)-1,1)
VarsData!1			= "case"
VarsData&case			= ispan(0,dimsizes(case)-1,1)
VarsData!2			= "lat"
VarsData&lat			= fspan(-89.5,89.5,nlat)
VarsData!3			= "lon"
VarsData&lon			= fspan(0.5,359.5,nlon)

do icase=0,dimsizes(case)-1
ftmp				= addfile(dir(icase)+case(icase)+append,"r")

do ivar=0,dimsizes(vars)-1
if(isfilevar(ftmp,vars(ivar)))then
	VarsData(ivar,icase,:,:)		= tofloat(ftmp->$vars(ivar)$(0,:,:))
end if

end do ; do ivar=

end do ; do icase=

printVarSummary(VarsData)

; weighting of latitude 

lat_here			= VarsData&lat
rad				= 4.0*atan(1.0)/180.
clat				= cos(lat_here*rad)
clat!0				= "lat"
clat&lat			= lat_here

lonS				= 0.
lonE				= 360.

; global mean
VarsData_GLB			= wgt_areaave_Wrap(VarsData(:,:,:,{lonS:lonE}),clat,1.0,1) 


print("------------------------------OBS--------------------------Park-RH------------------------Gauss-PDF---")
print(vars+"			"+sprintf("%9.3f",VarsData_GLB(:,0))+"		"+sprintf("%9.3f",VarsData_GLB(:,1))+"		"+sprintf("%9.3f",VarsData_GLB(:,2)))

exit
;=============================Make a table=================================
; Header
ncr1  = (/1,1/)               ; 1 row, 1 column
x1    = (/0.005,0.335/)       ; Start and end X
y1    = (/0.900,0.995/)       ; Start and end Y
text1 = "Metrics"

; Sub header
ncr2  = (/1,3/)               ; 1 rows, 18 columns
x2    = (/0.335,0.995/)       ; Start and end X
y2    = (/0.900,0.995/)       ; Start and end Y
text2 = (/"OBS","Park-RH","Gauss-PDF"/)

; cases index
ncr3	= (/dimsizes(vars),1/)
x3	= (/0.005,0.335/)
y3	= (/0.005,0.900/)

text3_tmp	= vars+" ("+vars_units+")"
text3		= transpose(text3_tmp)

; Main table body
ncr4  = (/dimsizes(vars),3/)              ; 5 rows, 17 columns
x4    = (/0.335,0.995/)       ; Start and end X
y4    = (/0.005,0.900/)       ; Start and end Y

text4_tmp		= (/ \
				sprintf("%9.3f",VarsData_GLB(:,0)),\
				sprintf("%9.3f",VarsData_GLB(:,1)),\
				sprintf("%9.3f",VarsData_GLB(:,2))\
			  /)
printVarSummary(text4_tmp)
text4				= transpose(text4_tmp)
printVarSummary(text4)


wks = gsn_open_wks("pdf","./figure/TOAbalance-Glbmean-table")  ; send graphics to PNG file

;
; Main header.
;
res1               = True
res1@txFontHeightF = 0.025
res1@gsFillColor   = "CornFlowerBlue"

gsn_table(wks,ncr1,x1,y1,text1,res1)

;
; Sub header
;
res2               = True
res2@txFontHeightF = 0.02
res2@gsFillColor   = "Gray"

gsn_table(wks,ncr2,x2,y2,text2,res2)

res3 = True             ; Set up resource list
res3@txFontHeightF = 0.02
res3@gsFillColor   = "Gray"

gsn_table(wks,ncr3,x3,y3,text3,res3)

;
; Main body of table.
;
res4 = True             ; Set up resource list
res4@txFontHeightF 	= 0.02

gsn_table(wks,ncr4,x4,y4,text4,res4)

frame(wks)      ; Advance the frame.

exit


