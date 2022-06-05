subroutine por(p,phi0,phi,phi_p)
    implicit none
    double precision, intent(in) :: p, phi0
    double precision, intent(out) :: phi, phi_p
    phi= phi0* exp(5.0*10.0**(-6.0)*(p-5000.0))
    phi_p= 5.0*10.0**(-6.0)*phi0* exp(5.0*10.0**(-6.0)*(p-5000.0))
end subroutine