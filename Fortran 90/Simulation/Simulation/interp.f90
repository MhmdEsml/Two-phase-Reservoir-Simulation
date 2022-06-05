subroutine interp(xx,yy,length,x,y)
    ! linear interpolation
    implicit none
    integer :: length
    double precision, dimension(1:length):: xx, yy
    double precision :: x, y
    integer :: i
        if (x < xx(1)) then
            y=yy(1)-(yy(1)-yy(2))*(xx(1)-x)/(xx(1)-xx(2))
            
        else if (x > xx(size(xx))) then
            y=yy(size(xx)-1)-(yy(size(xx)-1)-yy(size(xx)))*(xx(size(xx)-1)-x)/(xx(size(xx)-1)-xx(size(xx)))
        else
            do i=2,size(xx)
                if (x < xx(i)) then
                    y=yy(i-1)-(yy(i-1)-yy(i))*(xx(i-1)-x)/(xx(i-1)-xx(i))
                    exit
                end if
            end do
        end if
end subroutine
    