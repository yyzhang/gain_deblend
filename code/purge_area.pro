pro purge_area, ind, image, area, count_area=count_area

      count_maxima_ind=dblarr(n_elements(ind))+1.D
      image=double(image)
      ind=long(ind)
      n_ind=long(n_elements(ind))
      size_image=size(image)
      sizex=long(size_image[1])
      sizey=long(size_image[2])
      res=lonarr(sizex, sizey)-1000L
      ;print, 1
      test=call_external('./srcs/flood_shed.so','main', image, res, ind, n_ind, count_maxima_ind, sizex, sizey,$
                          verbose=verbose, show_all_output=verbose, /unload)
      ;print, 2
      ind_fil=where(count_maxima_ind ge  area)                 
      ind=ind[ind_fil]
      if keyword_set(count_area) then begin
         count_area=count_maxima_ind[ind_fil]
      endif
end