function flood_shed, imag, maxima_ind, count_maxima_ind
    size_image=size(imag)
    sizex=size_image[1]
    sizey=size_image[2]
    
    x_ind=dblarr(sizex, sizey)
    for i=0, sizex-1 do x_ind[i, *]=i
    y_ind=dblarr(sizex, sizey)
    for i=0, sizey-1 do y_ind[*, i]=i
    ind=where(x_ind ne 0 and x_ind ne sizex-1 and y_ind ne 0 and y_ind ne sizey-1)
    
    res=dblarr(sizex, sizey)-1000.0
    res[maxima_ind]=dindgen(n_elements(maxima_ind))
    
    
    ind_res=where(res lt -1)
    ind_res=ind_res(reverse(sort(imag[ind_res])))
    for i=0, n_elements(ind_res)-1 do begin
        ;print, i, n_elements(ind_res)
        x_ind_ij=x_ind[ind_res[i]]
        y_ind_ij=y_ind[ind_res[i]]
        ngb_xind=[x_ind_ij-1, x_ind_ij-1, x_ind_ij-1, x_ind_ij,   x_ind_ij,   x_ind_ij+1, x_ind_ij+1, x_ind_ij+1]
        ngb_yind=[y_ind_ij-1, y_ind_ij,   y_ind_ij+1, y_ind_ij-1, y_ind_ij+1, y_ind_ij-1, y_ind_ij,   y_ind_ij+1]
        ind_ngb=where(res[ngb_xind, ngb_yind] ge 0, count)
        if count gt 0 then begin
           idx=sort((res[ngb_xind, ngb_yind])(ind_ngb))
           indngb_uniq=uniq((res[ngb_xind, ngb_yind])(ind_ngb[idx]))
           if n_elements(indngb_uniq) eq 1 then begin
              res[ind_res[i]]=(res[ngb_xind, ngb_yind])(ind_ngb[indngb_uniq])
              count_maxima_ind[res[ind_res[i]]]=count_maxima_ind[res[ind_res[i]]]+1.0
           endif else res[ind_res[i]]=-1
        endif
    endfor
    
    
    return, res
end