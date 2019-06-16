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
season			=(/"ANN"/)
var3D			=(/"TEMP","UVEL","VVEL","WVEL"/)
var2D			=(/"TAUX","TAUY"/)
compare_type 		= (/"cmip-obs"/)
case			= (/"def","new"/)

Dir			= "/home/share3/lyl/work3/qinyi/double-ITCZ/data/"

ncdf = addfile("vars_2D_atm_cmip-obs.nc.nc" ,"r")
TAUX			= ncdf->TAUX
TAUY			= ncdf->TAUY

printVarSummary(TAUX)
printMinMax(TAUX,False)

;========================================================
; Prepare plot settings
;========================================================

latS			= -10.
latE			= -5.

lonS			= 160.
lonE			= 280.

LonValuePlot		= True


;========================================================
; Get the difference between two cases
;========================================================

vars_2D			= (/TAUX,TAUY/)
vars_2D(0,:,:,:)	= TAUX
vars_2D!0		= "vars"
vars_2D&vars		= (/0,1/)

printVarSummary(vars_2D)

vars_2D_diff		= vars_2D
vars_2D_diff(:,0,:,:)	= vars_2D(:,1,:,:)-vars_2D(:,0,:,:)
vars_2D_diff(:,1,:,:)	= vars_2D(:,2,:,:)-vars_2D(:,0,:,:)
vars_2D_diff(:,2,:,:)	= vars_2D(:,2,:,:)-vars_2D(:,1,:,:)
printVarSummary(vars_2D_diff)

;========================================================
; Get the zonal mean along zonal band
;========================================================

vars_2D_zm		= dim_avg_n_Wrap(vars_2D(:,:,{latS:latE},:),2)
printVarSummary(vars_2D_zm)

vars_2D_diff_zm		= dim_avg_n_Wrap(vars_2D_diff(:,:,{latS:latE},:),2)
printVarSummary(vars_2D_diff_zm)
printMinMax(vars_2D_diff_zm,True)


;---------------------------------------------------------------
; Draw Lon-Value Plots
;---------------------------------------------------------------
if(LonValuePlot)then
;---------------------------------------------------------------
; Draw Lon-Value Plots: separate schemes
;---------------------------------------------------------------
if(latE.lt.0)then
wks = gsn_open_wks("x11","./figure/ocn-heat-budget-LonValue-"+abs(latS)+"S-"+abs(latE)+"S")
else
wks = gsn_open_wks("x11","./figure/ocn-heat-budget-LonValue-"+abs(latS)+"S-"+abs(latE)+"N")
end if


plot				=new(2,graphic)

res 				= True

res@gsnFrame		 	= False
res@gsnDraw			= False

res@xyLineThicknessF 		= 3.0

res@trYReverse			= True

res@trYMinF			= 1.3*min(vars_2D_zm(0,:,{lonS:lonE}))
res@trYMaxF			= 0.7*max(vars_2D_zm(0,:,{lonS:lonE}))


;res@trXMinF			= lonS
;res@trXMaxF			= lonE

res@gsnLeftString		= "wind stress"
if(latE.lt.0)then
res@gsnRightString		= abs(latS)+"S - "+abs(latE)+"S"
else
res@gsnRightString		= abs(latS)+"S - "+abs(latE)+"N"
end if


res@tmXMajorGrid                = True          ; implement x grid 
res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
res@tmYMajorGrid                = True          ; implement y grid
res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmYMajorGridLineDashPattern = 2             ; select short dash lines

xsize				= dimsizes(vars_2D_zm(0,0,{lonS:lonE}))
res@tmXBMode   		= "Explicit"
res@tmXBValues 		= fspan(0,xsize,11)
res@tmXBLabels 		= (/"160E","","180","","160W","","140W","","120W","","100W"/)

res@xyLineColors		= (/"grey","red","blue"/)
res@xyDashPatterns		= (/0,0,0/)
plot(0)			 	= gsn_csm_y(wks,vars_2D_zm(0,:,{lonS:lonE}),res)

;========start to add legends
gres 				= True
gres@YPosPercent 		= 95.    ; expressed as %, 0->100, sets position of top border of legend 
                          		 ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
gres@XPosPercent 		= 5      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

lineres 			= True
lineres@lgLineColors 		= res@xyLineColors
lineres@lgLineThicknesses 	= res@xyLineThicknessF
lineres@lgDashIndexes		= res@xyDashPatterns

lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line


textres 			= True
textres@lgLabels 		= (/"OBS","DEF","NEW"/)

plot(0)				= simple_legend(wks,plot(0),gres,lineres,textres)

draw(plot(0))
frame(wks)
;---------------------------------------------------------------
; Draw Lon-Value Plots: differences
;---------------------------------------------------------------

res@trYMinF			= 1.3*min(vars_2D_diff_zm(0,:,{lonS:lonE}))
res@trYMaxF			= 1.3*max(vars_2D_diff_zm(0,:,{lonS:lonE}))


plot(1)			 	= gsn_csm_y(wks,vars_2D_diff_zm(0,:,{lonS:lonE}),res)

textres@lgLabels 		= (/"DEF-OBS","NEW-OBS","NEW-DEF"/)
plot(1)				= simple_legend(wks,plot(1),gres,lineres,textres)

draw(plot(1))
frame(wks)

gsn_panel(wks,plot,(/1,2/),False)

delete(res)
delete(plot)
delete(wks)

end if ; LonValuePlot

                               
end
