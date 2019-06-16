begin

  wks  = gsn_open_wks("x11","./figure/maponly")              ; send graphics to PNG file
;  plot = gsn_csm_map_ce(wks,False)                  ; draw default
;****************************************************
; now remove small islands from map
;****************************************************
  res                           = True              ; plot mods desired

  res@gsnDraw			= False
  res@gsnFrame			= False

  res@mpAreaMaskingOn           = True              ; turn on area masking
  res@mpMaskAreaSpecifiers      = "SmallIslands"    ; don't draw this area
  res@mpOutlineOn               = False             ; turn off outline

  res@mpMinLonF              = 80
  res@mpMaxLonF              = 300
  res@mpCenterLonF           = 180

  res@mpMinLatF			= -40.
  res@mpMaxLatF			= 0.

   ;res@tmXTOn = False ; Turn off tickmarks and labels 
   ;res@tmYROn = False 
   res@tmXBMinorOn	= False
   res@tmXTMinorOn	= False
   res@tmYLMinorOn	= False
   res@tmYRMinorOn	= False

   res@tmXBMajorLengthF	= 0.005
   res@tmXTMajorLengthF	= 0.005
   res@tmYLMajorLengthF	= 0.005
   res@tmYRMajorLengthF	= 0.005

   res@tmXBLabelFontHeightF	= 0.015
   res@tmYLLabelFontHeightF	= 0.015

   res@tmBorderThicknessF = 1.0

  plot = gsn_csm_map(wks,res)                    ; draw global map

;************************************************
; create points for box
;************************************************
 ypts = (/ -5.0,  -5.0,  -20.0,   -20.0, -5.0/)
 xpts = (/180.0, 230.0, 230.0, 180.0, 180.0/)
;************************************************
; add the box
;************************************************
  resp                  = True                      ; polyline mods desired
;  resp@gsLineColor      = "blue"                     ; color of lines
;  resp@gsLineThicknessF = 2.0                       ; thickness of lines
;  resp@gsLineLabelString= "test"                    ; adds a line label string

resp@gsFillColor	= "blue"
resp@gsFillOpacityF	= 0.2

  gon_1 = gsn_add_polygon(wks, plot, xpts, ypts, resp)

 ypts2 = (/ -5.0,  -5.0,  -20.0,   -20.0, -5.0/)
 xpts2 = (/230.0, 260.0, 260.0, 230.0, 230.0/)

 resp@gsFillColor	= "green"
 
 gon_2 = gsn_add_polygon(wks,plot,xpts2,ypts2,resp)


  draw(plot)
  frame(wks)


  end
