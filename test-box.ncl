begin

latSs				= (/-20,-20,-20,-20,-20,-20,-5,-5,-5/)
latEs				= (/-10,-10,-10,-5,-5,-5,5,5,5/)
lonSs				= (/180,230,260,180,230,260,240,200,160/)
lonEs				= (/225,255,285,225,255,285,275,235,195/)

latS				= -30.
latE				= 30.
lonS				= 150.
lonE				= 300.

nsteps				= 256/dimsizes(latSs)
colors				= ispan(2,256,nsteps)
dashes				= fspan(0,dimsizes(latEs)-1,dimsizes(latEs))


wks = gsn_open_wks("pdf","./figure/map-with-all-boxes")                 ; send graphics to PNG file

res           = True                              ; plot mods desired
res@gsnDraw   = False                             ; don't draw yet
res@gsnFrame  = False                             ; don't advance yet

res@mpMinLatF 				= latS
res@mpMaxLatF 				= latE
res@mpMinLonF 				= lonS
res@mpMaxLonF 				= lonE
res@mpCenterLonF			= 180.

  plot = gsn_csm_map(wks,res)          ; create the map plot

  res1					=	True
  res1@gsLineColor			=	"blue"
  res1@gsLineThicknessF			=	2.0
  res1@gsLineDashPattern		=	1

  equator	=	gsn_add_polyline(wks,plot,(/lonS,lonE/),(/0,0/),res1)


do iarea=0,dimsizes(latSs)-1
 ypts = (/ latSs(iarea),  latSs(iarea),  latEs(iarea),latEs(iarea), latSs(iarea)/)
 xpts = (/lonSs(iarea), lonEs(iarea),lonEs(iarea), lonSs(iarea),lonSs(iarea)/)


;************************************************
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = colors(iarea)
  resp@gsLineThicknessF = 3.0                       ; thickness of lines
  resp@gsLineDashPattern= dashes(iarea)

; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.

  dum = new(4,graphic)  

; draw each line separately. Each line must contain two points.
  do i = 0 , 3
    dum(i)=gsn_add_polyline(wks,plot,xpts(i:i+1),ypts(i:i+1),resp)      
  end do 

 str1 = unique_string("dum") 
 plot@$str1$ = dum

  draw(plot)
  delete(dum)
end do

  frame(wks)

end
