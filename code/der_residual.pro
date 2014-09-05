;+
; NAME:
; der_residual.pro
; 
; PURPOSE:
; Derive the light that are from blended astronomical sources by subtracting light from "background sources".
; Documentation can be found as Zhang et al., 2014, in preparation.
; 
; CALLING:
; res_str=der_residual(image, weight, maxima_ind, count_area)
; 
; INPUTS:
; image = 2-D array, weight=2-D array, maxima_ind, 1-d array of indices in image, 
; count_array=1-d array with same number of elements in maxima_ind.
; image is a astronomical image
; weight is the weight of this astronomical image.
; maxima_ind is a 1-d array of indices in image that are identified as blended astronomical sources. 
; count_array is the number of elements associated with each source.
; 
; KEYWORDS:
; 
; OUTPUTS:
; Function returns structure array with the same number of elements as maxima_ind and count_array.
; The structure has tags org_img, weight_img, subbkg_img, xcor_org, ycor_org, inpaint_radius and subbkg_flag.
; org_imag is a 51*51 pixels cut out of the orginal image for each blended source.
; wt_imag is a 51*51 pixels cut out of the orginal weight image for each blended source.
; subbkg_img is a 51*51 pixels stamp image of a blended source after subtracting light from "background sources".
; xcor_org and ycor_org is the x and y coordinates of a blended source in image.
; inpiant_radius is the inpainting radius used to subtract light.
; subbkg_flag: 1 "background sources" light subtraction performed. 0 "background sources" light subtraction NOT performed. 
; 
; HISTORY:
; Yuanyuan Zhang 2013
;-

function der_residual, image, weight, maxima_ind, seg_area
    half=25
    len=2*half+1
    str={org_img: dblarr(len, len), subbkg_img: dblarr(len, len), weight_img: dblarr(len, len), xcor_org:-1.D, ycor_org:-1.D, inpaint_radius: -1.D, subbkg_flag:0.D}
    str=replicate(str, n_elements(maxima_ind))
    if n_elements(maxima_ind) ne n_elements(seg_area) then print, 'Number of elements in argument[1] should equal that of argument[2]!' else begin
        size_image=size(image)
        sizex=long(size_image[1])
        sizey=long(size_image[2])
        x_ind=lonarr(sizex, sizey)
        for i=0, sizex-1 do x_ind[i, *]=i
        y_ind=lonarr(sizex, sizey)
        for i=0, sizey-1 do y_ind[*, i]=i
        
        if max(maxima_ind) gt -1 then begin
           for i=0, n_elements(maxima_ind)-1 do begin
               maxind=maxima_ind[i]
               ;print, i, n_elements(maxima_ind)
               if maxind gt -1 and maxind lt sizex*sizey then begin
                  xcor_org=x_ind(maxind)
                  ycor_org=y_ind(maxind)
                  if xcor_org gt half and xcor_org lt sizex-half and ycor_org gt half and ycor_org lt sizex-half then begin
                     org_img=double(image[xcor_org-half:xcor_org+half, ycor_org-half:ycor_org+half])
                     weight_img=double(weight[xcor_org-half:xcor_org+half, ycor_org-half:ycor_org+half])
                     grad_img=gradient(org_img, /vector)
                     gradx=org_img
                     gradx=double(smooth(grad_img[*, *, 0], 4))
                     grady=org_img
                     grady=double(smooth(grad_img[*, *, 1], 4))
                     grad_img=0
                     bkg_img=org_img
                     naxes1=long(len)
                     naxes2=long(len)
                     inpaint_radius=double(sqrt(seg_area[i]/!dpi)+1)
                     flag=1L
                     test=call_external('./srcs/fmm4.so','main', org_img, bkg_img, inpaint_radius, naxes1, naxes2, flag, gradx, grady, $
                              verbose=verbose, show_all_output=verbose, /unload)
                     str[i].org_img=org_img
                     str[i].subbkg_img=org_img-bkg_img
                     str[i].weight_img=weight_img
                     str[i].xcor_org=xcor_org
                     str[i].ycor_org=ycor_org
                     str[i].inpaint_radius=inpaint_radius
                     str[i].subbkg_flag=flag
                  endif else str[i].subbkg_flag=0
               endif else str[i].subbkg_flag=0
           endfor
        endif
        ;inf_filter=where(maxima_ind gt -1 and maxima_ind lt sizex*sizey)
        ;maxima_ind=maxima_ind(inf_filter)
        ;seg_area=seg_area(inf_filter)
        ;str=str[inf_filter]
    endelse
    return, str
end
