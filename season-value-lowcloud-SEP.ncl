

latS			= getenv("latS")
latE			= getenv("latE")
lonS			= getenv("lonS")
lonE			= getenv("lonE")

print(latS)

f1			= addfile("./CLDLOW_3cases_SEP_avg.nc","r")

data			= new((/3,12/),float)
data(0,:)		= f1->CLDLOW0
data(1,:)		= f1->CLDLOW1
data(2,:)		= f1->CLDLOW2

; change unit in model results
data(1:,:)		= data(1:,:)*100.

wks						= gsn_open_wks("pdf","./figure/CLDLOW_LonIs"+lonS+"_"+lonE+"_LatIs"+latS+"_"+latE)
res						= True

res@vpWidthF					= 0.7
res@vpHeightF					= 0.4

res@gsnFrame					= False
res@gsnDraw					= False

res@xyLineColors				= (/"black","red","blue"/)
res@xyDashPatterns				= (/0,0,0/)
res@xyLineThicknessF				= 3.

res@tiMainString				= "Lon="+lonS+" to "+lonE+"  Lat="+latS+" to "+latE

plot						= gsn_csm_y(wks,data,res)

legends                         		= (/"ISCCP","Park-RH","Gauss-PDF"/)

lgres                           	  	= True
lgres@DashIndexes				= res@xyDashPatterns
lgres@lgLineThicknessF 				= res@xyLineThicknessF
lgres@lgLineColors				= res@xyLineColors

lgres@lgLabelFontHeightF        		= 0.13
lgres@vpWidthF                  		= 0.20
lgres@vpHeightF                 		= 0.15
lgres@lgPerimOn                 		= False

lbid1                           		= gsn_create_legend(wks,dimsizes(legends),legends,lgres)

amres                           		= True
amres@amJust					= "TopLeft"
amres@amParallelPosF            		= -0.50
amres@amOrthogonalPosF          		= -0.50

annoid1                         		= gsn_add_annotation(plot,lbid1,amres)

draw(plot)
frame(wks)



