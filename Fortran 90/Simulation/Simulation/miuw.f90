subroutine miu_w (p,miuw)
    ! viscousity of water:
    double precision, intent(in) :: p
    double precision, intent(out) :: miuw
    double precision, dimension(1:2) :: visw, po
    visw=(/0.52341,0.56341/)
    po=(/3600.0,3900.0/)
    call interp(po, visw, size(po), p, miuw )
end subroutine