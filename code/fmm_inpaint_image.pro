function fmm_inpaint_image, image, impaint_f, epsilo ; impaint_f 1 area: to be inpainted. 0 area: known.
    imag=image
    size_image=size(image)
    impaint_t=dblarr(size_image[1], size_image[2])
    init_t=where(impaint_f gt 0)
    impaint_t[init_t]=10.0^6.0
    
    X_IND=DBLARR(SIZE_IMAGE[1], SIZE_IMAGE[2])
    FOR I=0, SIZE_IMAGE[1]-1 DO X_IND[I, *]=I
    Y_IND=DBLARR(SIZE_IMAGE[1], SIZE_IMAGE[2])
    FOR I=0, SIZE_IMAGE[2]-1 DO Y_IND[*, I]=I
    
    ;initialize impaint_f
    ind=where(impaint_f eq 0, count)
    ind=ind[(sort(impaint_t[ind]))]
    ind_ij=ind[0]
    while(count gt 0) do begin
         ind_ij=ind[count-1]
         impaint_f[ind_ij]=-1.0
         count=count-1
         
         x_ind_ij=x_ind[ind_ij]
         y_ind_ij=y_ind[ind_ij]
         ngb_xind=[x_ind_ij-1, x_ind_ij, x_ind_ij+1, x_ind_ij]
         ngb_yind=[y_ind_ij, y_ind_ij-1, y_ind_ij, y_ind_ij+1]
         ind_ngb=where(ngb_xind ge 0 and ngb_xind le size_image[1]-1 and ngb_yind ge 0 and ngb_yind le size_image[2]-1, count_ngb)
         if count_ngb gt 0 then begin
            ngb_xind=ngb_xind[ind_ngb]
            ngb_yind=ngb_yind[ind_ngb]
            for i=0, count_ngb-1 do begin
              if impaint_f[ngb_xind[i], ngb_yind[i]] eq 1 then impaint_f[ind_ij]=0.0
            endfor
         endif
    endwhile
    print, 'number of boudaries', n_elements(where(impaint_f lt 0.5 and impaint_f gt -0.5))
    
    Ttemp=image;temp_t(impaint_f, epsilo, image)
    impaint_t=ttemp;*alog(abs(tttemp))
    tttemp=image  
    
    ind=where(impaint_f eq 0, count)
    ;ind=ind[(sort(ttemp[ind]))]
    ind=ind[(sort(impaint_t[ind]))]
    ind_ij=ind[0]
    while(count gt 0) do begin
         impaint_f[ind_ij]=-1.0
         
         x_ind_ij=x_ind[ind_ij]
         y_ind_ij=y_ind[ind_ij]
         ngb_xind=[x_ind_ij-1, x_ind_ij, x_ind_ij+1, x_ind_ij]
         ngb_yind=[y_ind_ij, y_ind_ij-1, y_ind_ij, y_ind_ij+1]
         ind_ngb=where(ngb_xind ge 0 and ngb_xind le size_image[1]-1 and ngb_yind ge 0 and ngb_yind le size_image[2]-1, count_ngb)
         if count_ngb gt 0 then begin
            ngb_xind=ngb_xind[ind_ngb]
            ngb_yind=ngb_yind[ind_ngb]
            for i=0, count_ngb-1 do begin
              if impaint_f[ngb_xind[i], ngb_yind[i]] gt -0.5 then begin
                 if impaint_f[ngb_xind[i], ngb_yind[i]] eq 1 then begin; step 2
                    imag[ngb_xind[i], ngb_yind[i]]=impaint(ngb_xind[i], ngb_yind[i], imag, impaint_f, tttemp, epsilo)   ; step 3
                    impaint_f[ngb_xind[i], ngb_yind[i]]=0.0 
                 endif
                 ;impaint_t[ngb_xind[i], ngb_yind[i]]=;min([solve_t(impaint_t, impaint_f, ngb_xind[i]-1, ngb_yind[i], ngb_xind[i], ngb_yind[i]-1), $
                 ;                                         solve_t(impaint_t, impaint_f, ngb_xind[i]+1, ngb_yind[i], ngb_xind[i], ngb_yind[i]-1), $
                 ;                                         solve_t(impaint_t, impaint_f, ngb_xind[i]-1, ngb_yind[i], ngb_xind[i], ngb_yind[i]+1), $
                 ;                                         solve_t(impaint_t, impaint_f, ngb_xind[i]+1, ngb_yind[i], ngb_xind[i], ngb_yind[i]+1)])   ;step 4
              endif
            endfor
         endif
         ind=where(impaint_f eq 0, count) ; step 5 and step 1
         if count gt 0 then begin
            tttemp=imag
            ;print, impaint_t[ind[0]]
            ;ind=ind[(sort(ttemp[ind]))]
            ;mpaint_t=imag
            ind=ind[(sort(impaint_t[ind]))]
            ;print, impaint_t[ind[0]]
            ind_ij=ind[0]
         endif 
         ;print, ind_ij, count, impaint_f[ind_ij], n_elements(where(impaint_f eq -1)), n_elements(impaint_f)-n_elements(where(impaint_f eq -1))
    endwhile
    ;imag=exp(imag)-1.0+min_imag
    print, 'Greetings from fmm!'
    return, imag
    
    
end

function solve_t, t, f, i1, j1, i2, j2, f2
    sol=1.0*10.0^6.0
    size_image=size(t)
    size_1=size_image[1]
    size_2=size_image[2]
    if i1 ge 0 and i1 le size_1-1 and j1 ge 0 and j1 le size_2-1 then begin
    if f[i1, j1] eq -1 then begin
        ; if f(i1, j1)==known
        if i2 ge 0 and i2 le size_1-1 and j2 ge 0 and j2 le size_2-1 then begin
        if f[i2, j2] eq -1 then begin
           ; if f(i2, j2)==known
           r=2.0*f2-(t(i1, j1)-t(i2, j2))^2
           if r ge 0 then begin
              s=(t[i1, j1]+t[i2, j2]+sqrt(r))/2.0
              if s ge t[i1, j1] and s ge t[i2, j2] then sol=s 
           endif else sol=f2+t[i1, j1]
        endif else sol=f2+t[i1, j1]
        endif else sol=f2+t[i1, j1] ;if f(i2, j2) unknown, then extrapolate sol only from i1, j1
    endif else begin
        if i2 ge 0 and i2 le size_1-1 and j2 ge 0 and j2 le size_2-1 then begin
        if f[i2, j2] eq -1 then begin
           sol=f2+t[i2, j2] ;if only f(i2, j2) known, then extrapolate sol only from i2, j2
        endif
        endif
    endelse
    endif else begin
      if i2 ge 0 and i2 le size_1-1 and j2 ge 0 and j2 le size_2-1 then begin
        if f[i2, j2] eq -1 then begin
           sol=f2+t[i2, j2] ;if only f(i2, j2) known, then extrapolate sol only from i2, j2
        endif
      endif
    endelse
  
  return, sol
end

function impaint, i, j, image, f, t, epsilo
    ;imag=smooth(image, 3)
    Iij=-1;image[i, j]
    Ia=0.D
    s=0.D
    ;print, 'check image painted0:', Iij, image[i, j], image[i, j]-Iij
    size_image=size(image)
    x_ind=dblarr(size_image[1], size_image[2])
    for cs=0, size_image[1]-1 do x_ind[cs, *]=cs
    y_ind=dblarr(size_image[1], size_image[2])
    for cs=0, size_image[2]-1 do y_ind[*, cs]=cs
    distance=dblarr(size_image[1], size_image[2])
    distance=sqrt((x_ind-i)^2+(y_ind-j)^2)
    
    ;compute T around new boundary
    gradtij=grad(t, i, j)
    gradiij=grad(image, i, j)
    
    ind=where(distance le epsilo and f lt 0, count_ind)
    for counti=0, count_ind-1 do begin
        k=x_ind[ind[counti]]
        l=y_ind[ind[counti]]
        rvector=[i-k, j-l]
        
        r_len2=(total(rvector^2))
        if r_len2 gt 0.1 then begin
          
          dir=abs(total(rvector*gradtij))/sqrt(r_len2)
          dir=abs(total(rvector*gradiij))/sqrt(r_len2)
          dst=1.0/r_len2
          lev=1.0/(1.0+abs(t[k, l]-t[i,j]))
          ;lev=1.0/(1.0+abs(image[k, l]-image[i,j]))
          w=dst*lev;*dir;dst;*dir
          gradi=grad(image, k, l, f=f)
        
          Ia=Ia+w*(image[k, l]+total(gradi*rvector))
          s=s+w
          ;print,'grad', count, count_ind, gradi
          ;print, 'rvector', count, count_ind, rvector
          ;print, count, count_ind
          ;print, dir, dst, lev, w
          ;print, gradtij
          ;print, gradi
          ;print, rvector
          ;print, total(gradi*rvector)
          ;print, 'Ia, s, w', Ia, s, w, image[k, l], image[k, l]+total(gradi*rvector)
        endif
    endfor
    if s gt 0 then Iij=Ia/s
    ;print,'check image painted:' , Iij, image[i, j], image[i, j]-Iij 
    print, 'Ia, s, w', Ia, s, w
    if Iij lt min(image)-max(image) or Iij gt max(image)+max(image) then Iij=Iij
    return, Iij
end

function grad, T, i, j, f=f
   
   size_image=size(T)
   if keyword_set(f) eq 0 then begin
      f=dblarr(size_image[1], size_image[2])-1.0
   endif
   vecx=0.0
   vecy=0.0
   
   xlo=i-1
   xhi=i+1
   if xlo lt 0 then xlo=0
   if xhi gt size_image[1]-1 then xhi= size_image[1]-1
   ylo=j-1
   yhi=j+1
   if ylo lt 0 then ylo=0
   if yhi gt size_image[1]-1 then yhi= size_image[2]-1
   
   if f(xlo, j) gt 0.5 or f(xhi, j) gt 0.5 or f(i, ylo) gt 0.5 or f(i, yhi) gt 0.5 then begin
      xlo=i
      xhi=i
      ylo=j
      yhi=j
   endif
   
   vecx=T[xhi, j]-T[xlo, j]
   if abs(xhi-xlo) ge 0.01 then vecx=vecx/(xhi-xlo) else vecx=0.0
   vecy=T[i, yhi]-T[i, ylo]
   if abs(yhi-ylo) ge 0.01 then vecy=vecy/(yhi-ylo) else vecy=0.0
   vec=[vecx, vecy]
   
   return, vec
end

function temp_t, f, epsilo, image
    ; initialize impaint_t
    scale=max(image[where(f ge 0)])
    size_image=size(f)
    impaint_t=dblarr(size_image[1], size_image[2])
    X_IND=DBLARR(SIZE_IMAGE[1], SIZE_IMAGE[2])
    FOR I=0, SIZE_IMAGE[1]-1 DO X_IND[I, *]=I
    Y_IND=DBLARR(SIZE_IMAGE[1], SIZE_IMAGE[2])
    FOR I=0, SIZE_IMAGE[2]-1 DO Y_IND[*, I]=I
    
    ;calculate tout
    impaint_f=f
    tout=impaint_t
    ind_temp=where(f eq 1, count_temp)
    if count_temp gt 0 then impaint_f[ind_temp]=-1
    ind_temp=where(f eq -1, count_temp)
    if count_temp gt 0 then begin
       impaint_f[ind_temp]=1
       tout[ind_temp]=10.0^6.0
    endif
    
    ind=where(impaint_f eq 0, count)
    if count gt 0 then begin
       ind=ind[(sort(tout[ind]))]
       ind_ij=ind[0]
       min_t=tout[ind_ij]
    endif else min_t=10.0^6.0
    while(min_t le epsilo) do begin
         impaint_f[ind_ij]=-1.0
         
         x_ind_ij=x_ind[ind_ij]
         y_ind_ij=y_ind[ind_ij]
         ngb_xind=[x_ind_ij-1, x_ind_ij, x_ind_ij+1, x_ind_ij]
         ngb_yind=[y_ind_ij, y_ind_ij-1, y_ind_ij, y_ind_ij+1]
         ind_ngb=where(ngb_xind ge 0 and ngb_xind le size_image[1]-1 and ngb_yind ge 0 and ngb_yind le size_image[2]-1, count_ngb)
         if count_ngb gt 0 then begin
            ngb_xind=ngb_xind[ind_ngb]
            ngb_yind=ngb_yind[ind_ngb]
            for i=0, count_ngb-1 do begin
              if impaint_f[ngb_xind[i], ngb_yind[i]] gt -0.5 then begin
                 if impaint_f[ngb_xind[i], ngb_yind[i]] eq 1 then begin
                    impaint_f[ngb_xind[i], ngb_yind[i]]=0.0 ; step 2
                 endif
                 fij2=1;abs(image[ngb_xind[i], ngb_yind[i]]/scale)^2
                 tout[ngb_xind[i], ngb_yind[i]]=min([solve_t(tout, impaint_f, ngb_xind[i]-1, ngb_yind[i], ngb_xind[i], ngb_yind[i]-1, fij2), $
                                                          solve_t(tout, impaint_f, ngb_xind[i]+1, ngb_yind[i], ngb_xind[i], ngb_yind[i]-1, fij2), $
                                                          solve_t(tout, impaint_f, ngb_xind[i]-1, ngb_yind[i], ngb_xind[i], ngb_yind[i]+1, fij2), $
                                                          solve_t(tout, impaint_f, ngb_xind[i]+1, ngb_yind[i], ngb_xind[i], ngb_yind[i]+1, fij2)])   ;step 4
              endif
            endfor
         endif
         ind=where(impaint_f eq 0, count) ; step 5 and step 1
         if count gt 0 then begin
            ind=ind[(sort(tout[ind]))]
            ind_ij=ind[0]
            min_t=tout[ind_ij]
         endif else min_t=10.0^6.0
    endwhile
    tout=-tout
    
    ;calculate tin
    impaint_f=f
    tin=impaint_t
    ind_temp=where(f eq 1, count_temp)
    if count_temp gt 0 then begin
       tin[ind_temp]=10.0^6.0
    endif
    
    ind=where(impaint_f eq 0, count)
    if count gt 0 then begin
       ind=ind[(sort(tin[ind]))]
       ind_ij=ind[0]
    endif
    while(count gt 0) do begin
         impaint_f[ind_ij]=-1.0
         
         x_ind_ij=x_ind[ind_ij]
         y_ind_ij=y_ind[ind_ij]
         ngb_xind=[x_ind_ij-1, x_ind_ij, x_ind_ij+1, x_ind_ij]
         ngb_yind=[y_ind_ij, y_ind_ij-1, y_ind_ij, y_ind_ij+1]
         ind_ngb=where(ngb_xind ge 0 and ngb_xind le size_image[1]-1 and ngb_yind ge 0 and ngb_yind le size_image[2]-1, count_ngb)
         if count_ngb gt 0 then begin
            ngb_xind=ngb_xind[ind_ngb]
            ngb_yind=ngb_yind[ind_ngb]
            for i=0, count_ngb-1 do begin
              if impaint_f[ngb_xind[i], ngb_yind[i]] gt -0.5 then begin
                 if impaint_f[ngb_xind[i], ngb_yind[i]] eq 1 then begin
                    impaint_f[ngb_xind[i], ngb_yind[i]]=0.0 ; step 2
                 endif
                 fij2=1;abs(image[ngb_xind[i], ngb_yind[i]]/scale)^2
                 tin[ngb_xind[i], ngb_yind[i]]=min([solve_t(tin, impaint_f, ngb_xind[i]-1, ngb_yind[i], ngb_xind[i], ngb_yind[i]-1, fij2), $
                                                          solve_t(tin, impaint_f, ngb_xind[i]+1, ngb_yind[i], ngb_xind[i], ngb_yind[i]-1, fij2), $
                                                          solve_t(tin, impaint_f, ngb_xind[i]-1, ngb_yind[i], ngb_xind[i], ngb_yind[i]+1, fij2), $
                                                          solve_t(tin, impaint_f, ngb_xind[i]+1, ngb_yind[i], ngb_xind[i], ngb_yind[i]+1, fij2)])   ;step 4
              endif
            endfor
         endif
         ind=where(impaint_f eq 0, count)
         if count gt 0 then begin
            ind=ind[(sort(tin[ind]))]
            ind_ij=ind[0]
         endif
    endwhile
    
    impaint_t=tin
    ind=where(f eq 0, count)
    if count gt 0 then impaint_t[ind]=0.0
    ind=where(f lt -0.5, count)
    if count gt 0 then impaint_t[ind]=tout[ind]
    
    return, impaint_t
end