pro arrange_stamps1, res_str, temp_name, hdr_img, hdr_wt
    ma_ele=long(sqrt(n_elements(res_str))+1)
    size_ma=size(res_str.org_img)
    size_ma=long(size_ma[1])
    matrix0=dblarr(ma_ele*size_ma, ma_ele*size_ma)
    matrix1=dblarr(ma_ele*size_ma, ma_ele*size_ma)
    matrix2=dblarr(ma_ele*size_ma, ma_ele*size_ma)
    matrix3=dblarr(ma_ele*size_ma, ma_ele*size_ma)
    str={ind:-1L, xp:-1.0, yp:-1.0, xcor_org:-1.D, ycor_org:-1.D, inpaint_radius: -1.D, subbkg_flag:0.D}
    str=replicate(str, n_elements(res_str))
    for i=0, ma_ele-1 do begin
        for j=0, ma_ele-1 do begin
            ij=i+(j)*ma_ele
            if ij lt n_elements(res_str) then begin
               matrix0[i*size_ma:(i+1L)*size_ma-1, j*size_ma:(j+1L)*size_ma-1]=res_str[ij].subbkg_img
               matrix1[i*size_ma:(i+1L)*size_ma-1, j*size_ma:(j+1L)*size_ma-1]=res_str[ij].weight_img
               matrix2[i*size_ma:(i+1L)*size_ma-1, j*size_ma:(j+1L)*size_ma-1]=res_str[ij].org_img
               matrix3[i*size_ma:(i+1L)*size_ma-1, j*size_ma:(j+1L)*size_ma-1]=res_str[ij].org_img-res_str[ij].subbkg_img
               str[ij].ind=ij
               str[ij].xp=i*size_ma+(size_ma-1)/2L
               str[ij].yp=j*size_ma+(size_ma-1)/2L
               str[ij].xcor_org=res_str[ij].xcor_org
               str[ij].ycor_org=res_str[ij].ycor_org
               str[ij].inpaint_radius=res_str[ij].inpaint_radius
               str[ij].subbkg_flag=res_str[ij].subbkg_flag
            endif
        endfor
    endfor
    mwrfits, matrix0, temp_name, /create
    mwrfits, matrix1, temp_name
    mwrfits, matrix2, temp_name
    mwrfits, str, temp_name
    mwrfits, matrix3, temp_name
end