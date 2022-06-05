subroutine b_o_p (p,bo_p)
    ! derrivative of oil formation factor:
    double precision, intent(in) :: p
    double precision, intent(out) :: bo_p
    bo_p = 0.9802 * 10e-5 * exp( 10e-5 * (p-3600))
end subroutine