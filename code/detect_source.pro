;+
; NAME:
; detect_source.pro
; 
; PURPOSE:
; Identify maximas that look like genuine astronomical sources in an image, 
; and output number of pixels associated with these sources.
; Documentation can be found as Zhang et al., 2014, in preparation.
; 
; CALLING:
; maxima_ind=detect_source(image, count_area)
; 
; INPUTS:
; image = 2-D array
; 
; KEYWORDS:
; 
; OUTPUTS:
; Function returns 1-d array of indices of image that are identified as an astronomical source. 
; count_array returns number of elements associated with each source.
; 
; HISTORY:
; Yuanyuan Zhang 2013
;-

function detect_source, image, count_area
    size_image=size(image)
    sizex=size_image[1]
    sizey=size_image[2]
    
    min_image=min(image)
    log_image=image
    ind=where(log_image gt mean(image))
    log_image[ind]=alog(log_image[ind]/mean(image))+mean(image) ; change to log scale at above average to suppress high brightness contribution
    delta=(laplacian(image))
    exp_delta=(GAUSS_SMOOTH(delta, 1)); get the laplacian map
    img=GAUSS_SMOOTH(log_image*exp_delta, 1) ; convolve with image to produce the detection map.
    
    maxima_ind=maximas(img, 8); detect on detection map
    maxima_ind=maxima_ind[where(exp_delta[maxima_ind] gt mean(exp_delta)+0.0*stdev(exp_delta))] ; purge according to laplacian map
    maxima_ind=maxima_ind[where(img[maxima_ind] gt mean(img)+0.0*stdev(img))] ; purge according to detection map
    purge_area, maxima_ind, img,27 ; purge according to pixel numbers on detection map.
    print, 'Finished detecting from laplacian map.'
    img=0
    exp_delta=0
    delta=0
    log_image=0

    maxima_int=maximas(image, 8) ; local intensity peaks;
    ;help, maxima_int, /str
    ;print, max(maxima_int), min(maxima_int)
    purge_area, maxima_int, image, 27 ; purge according to pixel numbers.
    print, 'Finished detecting from intensity map.'
    
    sel_dis=2.0
    match=dblarr(n_elements(maxima_int))-1
    x_ind=dblarr(sizex, sizey)
    for i=0, sizex-1 do x_ind[i, *]=i
    y_ind=dblarr(sizex, sizey)
    for i=0, sizey-1 do y_ind[*, i]=i
    for i=0, n_elements(maxima_int)-1 do begin
        indi=maxima_int[i]
        dis=sqrt((x_ind[indi]-x_ind[maxima_ind])^2+(y_ind[indi]-y_ind[maxima_ind])^2)
        min_dis=min(dis)
        if min_dis lt sel_dis then match[i]=1
    endfor
    x_ind=0
    y_ind=0
    print, 'Finished matching the above two detections.'  
    maxima_int=maxima_int[where(match eq 1)]
    count_area=dblarr(n_elements(maxima_int))+1
    purge_area, maxima_int, image, 27, count_area=count_area; purge according to connected pixels above background
    
    print, 'Finished detecting!'    
    return, maxima_int
end
