;+
; NAME:
; gradient
; PURPOSE:
; Compute the gradient vector at every pixel of an image using
; neighbor pixels. Default is to return the Euclidean norm of gradient
; (2-D array), or specify keywords for other options.
; CALLING:
; gradim = gradient( image )
; INPUTS:
; image = 2-D array.
; KEYWORDS:
; XRANGE=, YRANGE= [min,max] range of pixel coordinates (for x,y axes),
;     allowing scaling of gradient magnitudes.
; /VECTOR then 2 images are returned (3-D array) with d/dx and d/dy.
; /ONE_SIDED then x & y-partial derivatives computed over adjacent pixels.
; /ABSVAL then magnitude of gradient is computed by absolute value norm.
; OUTPUTS:
; Function returns gradient norm (2-D array) or gradient vector (3-D).
; HISTORY:
; Frank Varosi U.Md. 1987
; F.V. NASA/GSFC 1991, added keyword options.
;-

function gradient, image, VECTOR=vector, ABSVAL_NORM=absval,  $
      ONE_SIDED=one_sided, XRANGE=xran, YRANGE=yran

  sim = size( image )

  if (sim(0) NE 2) then begin
    message,"expecting an image (2-D matrix)",/INFO
    return,sim
     endif

  if keyword_set( one_sided ) then begin
    didx = shift( image, -1,  0 ) - image
    didy = shift( image,  0, -1 ) - image
    endif else begin
    didx = ( shift( image, -1,  0 ) - shift( image, 1, 0 ) ) * 0.5
    didx(0,*) = image(1,*) - image(0,*)
    didy = ( shift( image,  0, -1 ) - shift( image, 0, 1 ) ) * 0.5
    didy(*,0) = image(*,1) - image(*,0)
     endelse

  didx(sim(1)-1,*) = image(sim(1)-1,*) - image(sim(1)-2,*)
  didy(*,sim(2)-1) = image(*,sim(2)-1) - image(*,sim(2)-2)

  if N_elements( xran ) EQ 2 then begin
    scale = sim(1) / float( xran(1)-xran(0) )
    didx = didx * scale
     endif

  if N_elements( yran ) EQ 2 then begin
    scale = sim(2) / float( yran(1)-yran(0) )
    didy = didy * scale
     endif

  if keyword_set( vector ) then  return, [ [[didx]], [[didy]] ]
  if keyword_set( absval ) then  return, abs( didx ) + abs( didy )

return, sqrt( didx*didx + didy*didy )
end