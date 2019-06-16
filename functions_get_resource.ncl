
function get_res_timelat(lat_min:numeric,lat_max:numeric)

begin

 res                  = True                ; plot mods desired
 res@gsnFrame = False
 res@gsnDraw = False
 res@cnFillOn         = True                ; color on
 res@lbOrientation    = "Vertical"          ; vertical label bar
 res@pmLabelBarOrthogonalPosF = -0.04       ; move label bar closer to axis
 res@cnLinesOn        = False               ; turn off contour lines

  res@tmXMajorGrid                = True          ; implement x grid 
  res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines

  res@tmYMinorGrid                = True          ; implement y grid
  res@tmYMinorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmYMinorGridLineDashPattern = 2             ; select short dash lines

 res@vpXF             = 0.12                ; default is 0.2 (aspect ratio)
 res@vpYF             = 0.8                 ; default is 0.8
 res@vpHeightF        = 0.4                 ; default is 0.6
 res@vpWidthF         = 0.75                ; default is 0.6

 res@tmXBMode = "Explicit"
 res@tmXBValues = fspan(0,11,12)
 res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

 res@lbLabelStride    = 2
 res@gsnStringFontHeightF = 0.025


return(res)

end



function get_res_latlon(lat_min:numeric,lat_max:numeric,lon_min:numeric,lon_max:numeric)
begin

res = True
res@gsnFrame = False
res@gsnDraw = False
res@gsnSpreadColors = True
res@cnLinesOn = False
res@cnFillOn = True

res@mpGridAndLimbOn        = True               ; turn on grid lines
res@mpGridLineDashPattern = 5
res@mpGridLatSpacingF = 10.
res@mpGridLonSpacingF = 20.
res@mpGridMaskMode = "MaskNotOcean"


;res@mpMinLonF = 0
;res@mpMaxLonF = 360
res@mpCenterLonF = 180

res@mpMinLatF =lat_min 
res@mpMaxLatF =lat_max

res@mpMinLonF = lon_min
res@mpMaxLonF = lon_max

res@tmXBLabelFontHeightF  = 0.014
res@tmYLLabelFontHeightF  = 0.014
; change label spacing to avoid overlap  

res@lbLabelStride    = 2                        ; every other label
res@tmXBTickSpacingF = 60.                      ; every 60 degrees 

res@lbOrientation   = "Vertical" 

res@gsnStringFontHeightF = 0.025

return(res)

end

function get_res_zonal_mean(lat_min:numeric,lat_max:numeric)
begin

res = True
res@gsnFrame = False
res@gsnDraw = False

res@trXMinF = lat_min
res@trXMaxF = lat_max

res@tmXBLabelFontHeightF  = 0.014
res@tmYLLabelFontHeightF  = 0.014

res@vpWidthF = 0.8
res@vpHeightF = 0.5

res@tmXMajorGrid                = True          ; implement x grid 
res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmXMajorGridLineDashPattern = 2             ; select short dash lines

res@tmYMajorGrid                = True          ; implement y grid
res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmYMajorGridLineDashPattern = 2             ; select short dash lines

;res@tmYMinorGrid                = True          ; implement y grid
;res@tmYMinorGridThicknessF      = 1.0           ; 2.0 is default
;res@tmYMinorGridLineDashPattern = 2             ; select short dash lines

return(res)



end
