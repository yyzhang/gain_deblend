pro wrapper_noncat, tile

    image_name='./example/eg_image.fits'
    image=mrdfits(image_name, 0, hdr_img, /silent)
    weight=mrdfits(image_name, 1, hdr_wt, /silent)
    size_image=size(image)
    sizex=long(size_image[1])
    sizey=long(size_image[2])
    
    cat_name='./example/eg_cat.fits'
    cat=mrdfits(cat_name, 1, /silent)
    maxima_ind=detect_source(image, count_area)

    external_x=long(cat.x_image)
    external_y=long(cat.y_image)
    purge_match, image, maxima_ind, count_area, external_x, external_y, 10

    res_str=der_residual(image, weight, maxima_ind, count_area) ; derive residual image using detection image
    tempres_name='./example/temp_res2.fits'
    arrange_stamps1, res_str, tempres_name, hdr_img, hdr_wt
    

    filename=tempres_name
    filepsf='./example/eg_image_psfcat.psf'
    cat_out='./example/eg_noncatdeb.fits'
    saturate_level=fxpar(hdr_img,'SATURATE')
    zeropoint=fxpar(hdr_img,'SEXMGZPT')

    ; SE call
    secommand=' sex '+filename+'[0],'+filename+'[0]'+$
                       ' -c ./example/sex.config '+$
                       ' -FILTER_NAME ./example/sex.conv '+$
                       ' -STARNNW_NAME ./example/sex.nnw '+$
                       ' -PARAMETERS_NAME ./example/sex.param_diskonly' +$
                       ' -CATALOG_TYPE FITS_1.0 '+$
                       ' -CATALOG_NAME '+cat_out+$
                       ' -WEIGHT_TYPE MAP_WEIGHT '+$
                       ' -WEIGHT_IMAGE '+filename+'[1],'+filename+'[1]'+$
                       ' -MEMOR_BUFFSIZE 2048 '+$
                       ' -MAG_ZEROPOINT '+strcompress(zeropoint,/remove_all)+$
                       ' -PSF_NAME '+filepsf+','+filepsf+$
                       ' -VERBOSE_TYPE NORMAL -DETECT_THRESH 1.5 -DEBLEND_MINCONT 0.005'
     spawn,secommand
        
     
     ; the following part insert the original posistion of de-blended sources into the SE catalog.  
     secat=mrdfits(cat_out, 1, /silent)
     str=secat[0]
     tcat2=mrdfits(tempres_name, 3, /silent)
     str=replicate(str, n_elements(tcat2))
     match_ind=lonarr(n_elements(tcat2))-1
     for ti=0, n_elements(tcat2)-1 do begin
         dis=sqrt((secat.x_image-tcat2[ti].xp)^2+(secat.y_image-tcat2[ti].yp)^2)
         min_dis=min(dis, min_ind)
         if min_dis lt 10 and min_ind gt -1 then begin
            str[ti]=secat[min_ind]
            match_ind[ti]=min_ind
         endif
      endfor
      str2=create_struct(str[0], 'sematch_flag', 0L, 'xp', -1.D, 'yp', -1.D, 'xcor_org', -1.D, 'ycor_org', -1.D, 'inpaint_radius', 0.D, 'subbkg_flag', 0.D, 'ind', -1L)
      str2=replicate(str2, n_elements(str))
      struct_assign, str, str2
      str=str2
      indcat2=where( match_ind gt -1)
      str[indcat2].sematch_flag=1
      str[indcat2].xp=tcat2[indcat2].xp
      str[indcat2].yp=tcat2[indcat2].yp
      str[indcat2].xcor_org=tcat2[indcat2].xcor_org
      str[indcat2].ycor_org=tcat2[indcat2].ycor_org
      str[indcat2].inpaint_radius=tcat2[indcat2].inpaint_radius
      str[indcat2].subbkg_flag=tcat2[indcat2].subbkg_flag
      str[indcat2].ind=tcat2[indcat2].ind
      mwrfits, str, cat_out, /create

end
