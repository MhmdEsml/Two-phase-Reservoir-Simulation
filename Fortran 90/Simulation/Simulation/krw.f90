subroutine kr_w (sw, krw, krw_s)
    ! relative perm of water:
    double precision, intent(in) :: sw
    double precision, intent(out) :: krw, krw_s
    ! swc=0.15, sor=0.3, krw_0=0.5, a=3
    ! krw= krw_0*(1-sw-swc)/(1-swc-sor)
    krw= 0.5* (sw-0.15)**3.0/(1.0-0.15-0.3)
    krw_s= -0.5*3* (sw-0.15)**2.0/(1.0-0.15-0.3)
    if (krw<0) then
        krw=0
        krw_s=0
    end if
end subroutine