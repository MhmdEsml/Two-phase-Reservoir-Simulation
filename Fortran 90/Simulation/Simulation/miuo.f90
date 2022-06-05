subroutine miu_o (p,miuo)
    ! viscousity of water:
    double precision, intent(in) :: p
    double precision, intent(out) :: miuo
    double precision, dimension(1:8) :: viso, po
    viso=(/1.16, 1.164, 1.167, 1.172, 1.177, 1.181, 1.185, 1.190/)
    po=(/400.0, 1200.0, 2000.0, 2800.0, 3600.0, 4400.0, 5200.0, 5600.0/)
    call interp(po, viso, size(po), p, miuo )
end subroutine