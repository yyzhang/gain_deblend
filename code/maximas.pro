function maximas, imag, area
    size_image=size(imag)
    sizex=long(size_image[1])
    sizey=long(size_image[2])
    
    res=lonarr(sizex, sizey)-1
    area=long(area)
    imag=double(imag)
    
    test=call_external('./srcs/maxima_2d.so','main', imag, res, area, sizex, sizey,$
                       verbose=verbose, show_all_output=verbose, /unload)
    maxima_ind=where(res gt 0)
    return, maxima_ind
end