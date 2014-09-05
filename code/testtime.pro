pro testtime
    image_name='./example/eg_image.fits'
    cat_name='./example/eg_cat.fits'
    cat=mrdfits(cat_name, 1)
    image=mrdfits(image_name, 0, hdr_img, /silent)
    weight=mrdfits(image_name, 1, hdr_wt, /silent)

    external_x=cat.x_image
    external_y=cat.y_image
    
    
    
    ; this block detect for sources
    t0=systime(/seconds)
    maxima_ind=detect_source(image, count_area)
    t1=systime(/seconds)
    print, 'Detection takes ', strtrim(string(t1-t0, format='(I20)'), 2), ' seconds.'
    
    ;this block purgest sources that are already in user-supplied catalog. 
    purge_match, image, maxima_ind, count_area, external_x, external_y, 10
    t2=systime(/seconds)
    print, 'Purging sources that are already in catalogs takes ', strtrim(string(t2-t1, format='(I20)'), 2), ' seconds.'

    ;this block derives residulas for sources.
    res_str=der_residual(image, weight, maxima_ind, count_area)
    t3=systime(/seconds)
    print, 'Deriving residual for sources takes ', strtrim(string(t3-t2, format='(I20)'), 2), ' seconds.'

end
