;+
; NAME:
; purge_match.pro
; 
; PURPOSE:
; Match the GAIN source list to the user supplied source list, 
; and single out those sources not identified by the original reduction as a list of newly identified, deblended sources.
; Documentation can be found as Zhang et al., 2014, in preparation.
; 
; CALLING:
; purge_match, image, maxima_ind, count_area, external_x, external_y, match_dis
; 
; INPUTS:
; image = 2-D array, maxima_ind=1-D array, count_area=1-D array, external_x= 1-D array, external_y=1-d array, match_dis=double
; maxima_ind is a 1-d array of indices in image that are identified as astronomical sources. 
; count_array is the number of elements associated with each source.
; external_x and external_y are lists of user-supplied x, y coordinates of image. 
; The astronomical sources at (external_x, external_y) will be purged out from maxima_ind and count_area.
; match_dis is the maximum distance allowed to match a source in maxima_ind to sources at (external_x, external_y). 
; 
; OUTPUTS:
; maxima_ind, count_area
; maxima_ind is a 1-d array of indices in image that are identified as astronomical sources after purging sources in the user-supplied list. 
; count_array is the number of elements associated with each source.
; 
; HISTORY:
; Yuanyuan Zhang 2013
;-

pro purge_match, image, maxima_ind, count_area, external_x, external_y, match_dis
    size_image=size(image)
    sizex=long(size_image[1])
    sizey=long(size_image[2])
    ind=where(external_x gt 0 and external_x lt sizex and external_y gt 0 and external_y lt sizey )
    external_x=long(external_x[ind])
    external_y=long(external_y[ind])
    external_indices=(external_x+external_y*sizex)
    
    external_indices=external_indices[reverse(sort(image[external_indices]))]
    external_indices=long(external_indices)
    maxima_ind=long(maxima_ind)
    image=double(image)
    nex=long(n_elements(external_indices))
    nma=long(n_elements(maxima_ind))
    match_dis=double(match_dis)
    sta=call_external('./srcs/purge_match.so','main', image, external_indices, nex, maxima_ind, nma, sizex, sizey, match_dis, $
                          verbose=verbose, show_all_output=verbose, /unload)
    ind=where(maxima_ind gt -0.5)
    maxima_ind=maxima_ind[ind]
    count_area=count_area[ind]
    
end