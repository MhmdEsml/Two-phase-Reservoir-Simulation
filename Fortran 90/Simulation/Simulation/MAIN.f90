program main
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! declration !!!!!!!!!!!!!!!!

    double precision :: alpha = 5.615, beta = 1.127, error, start, finish
    double precision :: miuo,miuo1, miuw, miuw1, kro, kro1, krw, krw1, kro_s, kro_s1, krw_s ,krw_s1 , bo, bon, inv_bo_p,&
        bo1, bw, bwn, inv_bw_p, bw1, phi, phi_p, to1, tw1, cop, cwp, cow, cww
    integer, parameter :: dt=5, nt=10, nx=15, ny=15, nz=1, dim=15*15*1
    double precision, parameter :: dx=75.0, dy=75.0, dz=30.0, vb=75.0*75.0*30.0, ax=75.0*30.0, ay=75.0*30.0, az=75.0*75.0
    double precision, parameter :: z0=-1000.0, gamma_o=0.7, gamma_w=1.0;
    double precision, dimension(1:dim) :: kx, ky, kz, phi0, z, pon, swn, po, sw, qw, qo
    double precision, dimension(1:2*dim,1) :: r, x
    double precision, dimension(1:2*dim,1:2*dim) :: jac
    integer :: i, j, k, t, m, n, p, INFO
    double precision, dimension(1,2*dim) :: IPIV
    
    qw(225)=10
    qo(1)=-10

    open(unit=1, file='por.txt', status='old', action='read' )
    read(1,*) phi0
    close(1)
    open(unit=2, file='perm_h.txt', status='old', action='read' )
    read(2,*) kx
    close(2)
    open(unit=3, file='perm_h.txt', status='old', action='read' )
    read(3,*) ky
    close(3)
    open(unit=4, file='perm_v.txt', status='old', action='read' )
    read(4,*) kz
    close(4)

    do k=1, nz
        do n=1, nx*ny
            z(n+nx*ny*(k-1))= z0-(k-1)*dz
        end do
    end do
    
    do i=1, dim
        pon(i)=5000.0
        swn(i)=0.3
    end do
    po=pon
    sw=swn
    
    call cpu_time(start)
    do t=1,nt
        error=1
        do while (error>1.0e-3)
            n=0
            
            do m=1, 2*dim
                r(m,1)=0
                do p=1, 2*dim
                    jac(m,p)=0
                end do
            end do
            
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        n=n+1
                        call miu_o(po(n), miuo)
                        call miu_w(po(n), miuw)
                        call kr_o(sw(n), kro, kro_s)
                        call kr_w(sw(n), krw, krw_s)
                        call b_o(pon(n), bon, inv_bo_p)
                        call b_o(po(n), bo, inv_bo_p)
                        call b_w(pon(n), bwn, inv_bw_p)
                        call b_w(po(n), bw, inv_bw_p)
                        call por(po(n),phi0(n),phi,phi_p)
                        
                        cop=(1-swn(n))*vb*(phi_p/bon+phi*inv_bo_p)/(alpha*dt)
                        cwp=swn(n)*vb*(phi_p/bwn+phi*inv_bw_p)/(alpha*dt)
                        cow=-vb*(phi/bo)/(alpha*dt)
                        cww=vb*(phi/bw)/(alpha*dt)
                        
                        r(n,1)=-cop*(po(n)-pon(n))-cow*(sw(n)-swn(n))+qo(n)
                        r(n+dim,1)=-cwp*(po(n)-pon(n))-cww*(sw(n)-swn(n))+qw(n)
                        
                        jac(n,n)=-cop
                        jac(n,n+dim)=-cow
                        
                        jac(n+dim,n)=-cwp
                        jac(n+dim,n+dim)=-cww
                        
                        ! -x:
                        if (i /= 1) then
                            call miu_o(po(n-1), miuo1)
                            call miu_w(po(n-1), miuw1)
                            call kr_o(sw(n-1), kro1, kro_s1)
                            call kr_w(sw(n-1), krw1, krw_s1)
                            call b_o(po(n-1), bo1, inv_bo_p)
                            call b_w(po(n-1), bw1, inv_bw_p)
                            if (po(n)>=po(n-1)) then
                                to1=beta*0.5*ax*(kx(n)+kx(n-1))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro
                                tw1=beta*0.5*ax*(kx(n)+kx(n-1))/dx* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw
                                
                                r(n,1)=r(n,1)+to1*(po(n-1)-po(n)-gamma_o*(z(n-1)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n-1)-po(n)-gamma_w*(z(n-1)-z(n)))
                             
                                jac(n, n)=jac(n,n)-to1 ! oil r/p(n)
                                jac(n,n-1)=to1 ! oil r/p(n-1)
                                
                                jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ax*(kx(n)+kx(n-1))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n-1)-po(n)-gamma_o*(z(n-1)-z(n))) ! oil r/sw(n)
                                                            
                                jac(n+dim, n)=jac(n+dim, n)-tw1 ! water r/p(n)
                                jac(n+dim, n-1)=tw1 ! water r/p(n-1)
                            
                                jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*ax*(kx(n)+kx(n-1))/dx* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n-1)-po(n)-gamma_w*(z(n-1)-z(n))) ! water r/sw(n)
                                
                            else
                                to1=beta*0.5*ax*(kx(n)+kx(n-1))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1
                                tw1=beta*0.5*ax*(kx(n)+kx(n-1))/dx* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1
                                
                                r(n,1)=r(n,1)+to1*(po(n-1)-po(n)-gamma_o*(z(n-1)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n-1)-po(n)-gamma_w*(z(n-1)-z(n)))
                                
                                jac(n, n)=jac(n,n)-to1 ! oil r/p(n)
                                jac(n,n-1)=to1 ! oil r/p(n-1)
                            
                                jac(n, n+dim-1)=beta*0.5*ax*(kx(n)+kx(n-1))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n-1)-po(n)-gamma_o*(z(n-1)-z(n))) ! oil r/sw(n-1)
                                
                                jac(n+dim, n)=jac(n+dim, n)-tw1 ! water r/p(n)
                                jac(n+dim, n-1)=tw1 ! water r/p(n-1)
                            
                                jac(n+dim, n+dim-1)=beta*0.5*ax*(kx(n)+kx(n-1))/dx* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n-1)-po(n)-gamma_w*(z(n-1)-z(n))) ! water r/sw(n-1)
                                
                           end if                            
                        end if
                        
                        ! +x:
                        if (i /= nx) then
                            call miu_o(po(n+1), miuo1)
                            call miu_w(po(n+1), miuw1)
                            call kr_o(sw(n+1), kro1, kro_s1)
                            call kr_w(sw(n+1), krw1, krw_s1)
                            call b_o(po(n+1), bo1, inv_bo_p)
                            call b_w(po(n+1), bw1, inv_bw_p)
                            
                            if (po(n)>=po(n+1)) then
                                to1=beta*0.5*ax*(kx(n)+kx(n+1))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro
                                tw1=beta*0.5*ax*(kx(n)+kx(n+1))/dx* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw
                            
                                r(n,1)=r(n,1)+to1*(po(n+1)-po(n)-gamma_o*(z(n+1)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n+1)-po(n)-gamma_w*(z(n+1)-z(n)))
                            
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n+1)=to1
                            
                                jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ax*(kx(n)+kx(n+1))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n+1)-po(n)-gamma_o*(z(n+1)-z(n)))
                                                            
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n+1)=tw1
                            
                                jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*ax*(kx(n)+kx(n+1))/dx* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n+1)-po(n)-gamma_w*(z(n+1)-z(n)))
                                
                            else
                                to1=beta*0.5*ax*(kx(n)+kx(n+1))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1
                                tw1=beta*0.5*ax*(kx(n)+kx(n+1))/dx* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1
                            
                                r(n,1)=r(n,1)+to1*(po(n+1)-po(n)-gamma_o*(z(n+1)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n+1)-po(n)-gamma_w*(z(n+1)-z(n)))
                            
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n+1)=to1
                            
                                jac(n, n+dim+1)=beta*0.5*ax*(kx(n)+kx(n+1))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n+1)-po(n)-gamma_o*(z(n+1)-z(n)))
                            
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n+1)=tw1
                            
                                jac(n+dim, n+dim+1)=beta*0.5*ax*(kx(n)+kx(n+1))/dx* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n+1)-po(n)-gamma_w*(z(n+1)-z(n)))
                            end if                                
                        end if
                        
                        ! -y:
                        if (j /= 1) then
                            call miu_o(po(n-nx), miuo1)
                            call miu_w(po(n-nx), miuw1)
                            call kr_o(sw(n-nx), kro1, kro_s1)
                            call kr_w(sw(n-nx), krw1, krw_s1)
                            call b_o(po(n-nx), bo1, inv_bo_p)
                            call b_w(po(n-nx), bw1, inv_bw_p)
                            if (po(n)>=po(n-nx)) then
                                to1=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro
                                tw1=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw
                                
                                r(n,1)=r(n,1)+to1*(po(n-nx)-po(n)-gamma_o*(z(n-nx)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n-nx)-po(n)-gamma_w*(z(n-nx)-z(n)))
                                
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n-nx)=to1
                                
                                jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ay*(ky(n)+ky(n-nx))/dy* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n-nx)-po(n)-gamma_o*(z(n-nx)-z(n)))
                                
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n-nx)=tw1
                                
                                jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*ay*(ky(n)+ky(n-nx))/dy* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n-nx)-po(n)-gamma_w*(z(n-nx)-z(n)))
                            else
                                to1=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1
                                tw1=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1
                                
                                r(n,1)=r(n,1)+to1*(po(n-nx)-po(n)-gamma_o*(z(n-nx)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n-nx)-po(n)-gamma_w*(z(n-nx)-z(n)))
                                
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n-nx)=to1
                                
                                jac(n, n+dim-nx)=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n-nx)-po(n)-gamma_o*(z(n-nx)-z(n)))
                                
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n-nx)=tw1
                                
                                jac(n+dim, n+dim-nx)=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* &
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n-nx)-po(n)-gamma_w*(z(n-nx)-z(n)))
                            end if
                        end if
                        
                        ! +y:
                        if (j /= ny) then
                            call miu_o(po(n+nx), miuo1)
                            call miu_w(po(n+nx), miuw1)
                            call kr_o(sw(n+nx), kro1, kro_s1)
                            call kr_w(sw(n+nx), krw1, krw_s1)
                            call b_o(po(n+nx), bo1, inv_bo_p)
                            call b_w(po(n+nx), bw1, inv_bw_p)
                            if (po(n)>=po(n+nx)) then
                                to1=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro

                                tw1=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw
                            
                                r(n,1)=r(n,1)+to1*(po(n+nx)-po(n)-gamma_o*(z(n+nx)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n+nx)-po(n)-gamma_w*(z(n+nx)-z(n)))
                            
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n+nx)=to1
                            
                                jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ay*(ky(n)+ky(n+nx))/dy* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n+nx)-po(n)-gamma_o*(z(n+nx)-z(n)))
                            
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n+nx)=tw1
                            
                                jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*ay*(ky(n)+ky(n+nx))/dy* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n+nx)-po(n)-gamma_w*(z(n+nx)-z(n)))
                                
                            else
                                to1=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1
                                tw1=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1
                            
                                r(n,1)=r(n,1)+to1*(po(n+nx)-po(n)-gamma_o*(z(n+nx)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n+nx)-po(n)-gamma_w*(z(n+nx)-z(n)))
                            
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n+nx)=to1
                            
                                jac(n, n+dim+nx)=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n+nx)-po(n)-gamma_o*(z(n+nx)-z(n)))
                            
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n+nx)=tw1
                            
                                jac(n+dim, n+dim+nx)=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n+nx)-po(n)-gamma_w*(z(n+nx)-z(n)))
                            end if   
                        end if
                        
                        ! -z:
                        if (k /= 1) then
                            call miu_o(po(n-nx*ny), miuo1)
                            call miu_w(po(n-nx*ny), miuw1)
                            call kr_o(sw(n-nx*ny), kro1, kro_s1)
                            call kr_w(sw(n-nx*ny), krw1, krw_s1)
                            call b_o(po(n-nx*ny), bo1, inv_bo_p)
                            call b_w(po(n-nx*ny), bw1, inv_bw_p)
                            if (po(n)>=po(n-nx*ny)) then
                                to1=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro
                                tw1=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw
                            
                                r(n,1)=r(n,1)+to1*(po(n-nx*ny)-po(n)-gamma_o*(z(n-nx*ny)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n-nx*ny)-po(n)-gamma_w*(z(n-nx*ny)-z(n)))
                            
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n-nx*ny)=to1
                            
                                jac(n, n+dim)=jac(n, n+dim)+beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n-nx*ny)-po(n)-gamma_o*(z(n-nx*ny)-z(n)))
                            
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n-nx*ny)=tw1
                            
                                jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n-nx*ny)-po(n)-gamma_w*(z(n-nx*ny)-z(n)))

                            else
                                to1=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1
                                tw1=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1
                            
                                r(n,1)=r(n,1)+to1*(po(n-nx*ny)-po(n)-gamma_o*(z(n-nx*ny)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n-nx*ny)-po(n)-gamma_w*(z(n-nx*ny)-z(n)))
                            
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n-nx*ny)=to1
                            
                                jac(n, n+dim-nx*ny)=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n-nx*ny)-po(n)-gamma_o*(z(n-nx*ny)-z(n)))
                            
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n-nx*ny)=tw1
                            
                                jac(n+dim, n+dim-nx*ny)=beta*0.5*az*(kz(n)+kz(n-nx))/dz* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n-nx*ny)-po(n)-gamma_w*(z(n-nx*ny)-z(n)))
                            end if
                            
                            
                        end if
                        
                        ! +z:
                        if (k /= nz) then
                            call miu_o(po(n+nx*ny), miuo1)
                            call miu_w(po(n+nx*ny), miuw1)
                            call kr_o(sw(n+nx*ny), kro1, kro_s1)
                            call kr_w(sw(n+nx*ny), krw1, krw_s1)
                            call b_o(po(n+nx*ny), bo1, inv_bo_p)
                            call b_w(po(n+nx*ny), bw1, inv_bw_p)
                            if (po(n)>=po(n+nx*ny)) then
                                to1=beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro
                                tw1=beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw
                            
                                r(n,1)=r(n,1)+to1*(po(n+nx*ny)-po(n)-gamma_o*(z(n+nx*ny)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n+nx*ny)-po(n)-gamma_w*(z(n+nx*ny)-z(n)))
                            
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n+nx*ny)=to1
                            
                                jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ax*(kx(n)+kx(n+nx*ny))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n+nx*ny)-po(n)-gamma_o*(z(n+nx*ny)-z(n)))
                            
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n+nx*ny)=tw1
                            
                                jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n+nx*ny)-po(n)-gamma_w*(z(n+nx*ny)-z(n)))

                            else
                                to1=beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1
                                tw1=beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1
                            
                                r(n,1)=r(n,1)+to1*(po(n+nx*ny)-po(n)-gamma_o*(z(n+nx*ny)-z(n)))
                                r(n+dim,1)=r(n+dim,1)+tw1*(po(n+nx*ny)-po(n)-gamma_w*(z(n+nx*ny)-z(n)))
                            
                                jac(n, n)=jac(n, n)-to1
                                jac(n,n+nx*ny)=to1
                            
                                jac(n, n+dim+nx*ny)=beta*0.5*ax*(kx(n)+kx(n+nx*ny))/dx* &
                                    0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n+nx*ny)-po(n)-gamma_o*(z(n+nx*ny)-z(n)))
                            
                                jac(n+dim, n)=jac(n+dim, n)-tw1
                                jac(n+dim, n+nx*ny)=tw1
                            
                                jac(n+dim, n+dim+nx*ny)=beta*0.5*az*(kz(n)+kz(n-nx))/dz* &
                                    0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n+nx*ny)-po(n)-gamma_w*(z(n+nx*ny)-z(n)))
                            end if
                            
                        end if
  
                    end do ! i
                end do ! j
            end do ! k
            
            do p=1, 2*dim
                r(p,1)=-r(p,1)
            end do

            CALL DGETRF( 2*dim, 2*dim , jac, 2*dim , IPIV, INFO ) 
            CALL DGETRS('N' ,  2*dim , 1 , jac, 2*dim, IPIV, r , 2*dim , INFO)

    
            do p=1,dim
                po(p)=po(p)+r(p,1)
                sw(p)=sw(p)+r(p+dim,1)
            end do
            
            error= dnrm2(450,r,1)
            ! print*, error
            !error=0.0001
        end do ! while
        pon=po
        swn=sw
    end do ! timesteps
    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start
    print*, sw,po
    pause
end ! main
    