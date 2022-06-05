subroutine b_o (p,bo,inv_bo_p)
    ! oil formation factor:
    double precision, intent(in) :: p
    double precision, intent(out) :: bo, inv_bo_p
    bo = 0.9802 * exp( -10.0**(-5.0) * (p-3600.0))
    inv_bo_p = (1/0.9802) * 10.0**(-5.0) * exp( 10.0**(-5.0) * (p-3600.0))
end subroutine