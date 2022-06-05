clear, clc

alpha = 5.615; beta = 1.127;
dt=5; nt=10; nx=15; ny=15; nz=1; dim=nx*ny*nz;
dx=75.0; dy=75.0; dz=30.0; vb=dx*dy*dz; ax=dy*dz; ay=dx*dz; az=ax*ay;
z0=-1000.0; gamma_o=0.7; gamma_w=1.0;

qo=zeros(1,dim); qw=zeros(1,dim);
qo(1)=-10;
% qw(113)=-1.0;
% qw(1)=-50;
qw(225)=10;

phi0=importdata('por.mat');
kx=importdata('perm_h.mat');
ky=importdata('perm_h.mat');
kz=importdata('perm_v.mat');

for k=1:nz
    for n=1:nx*ny
        z(n+nx*ny*(k-1))= z0-(k-1)*dz;
    end
end

for i=1:dim
    pon(i)=5000;
    swn(i)=0.3;
end

po=pon;
sw=swn;

time=cputime;
for t=1:nt
    error=1;
    while error>10^-3
        n=0;
        
        r=zeros(2*dim,1);
        jac=zeros(2*dim);
        
        for k=1:nz
            for j=1:ny
                for i=1:nx
                    n=n+1;
                    
                    miuo= miu_o(po(n));
                    miuw= miu_w(po(n));
                    [kro,kro_s]=kr_o(sw(n));
                    [krw,krw_s]=kr_w(sw(n));
                    [bon, ~]= b_o(pon(n));
                    [bo, inv_bo_p]=b_o(po(n));
                    [bwn, ~]=b_w(pon(n));
                    [bw, inv_bw_p]=b_w(po(n));
                    [phi,phi_p]=poro(po(n),phi0(n));
                    
                    cop=(1-swn(n))*vb*(phi_p/bon+phi*inv_bo_p)/(alpha*dt);
%                    cop=(1-swn(n))*vb*phi_p/(alpha*dt*bo); %
                    cwp=swn(n)*vb*(phi_p/bwn+phi*inv_bw_p)/(alpha*dt);
%                    cwp=swn(n)*vb*phi_p/(alpha*dt*bw); %
                    cow=-vb*(phi/bo)/(alpha*dt);
%                    cow=-vb*phi/(alpha*dt*bo); %
                    cww=vb*(phi/bw)/(alpha*dt);
%                    cww=vb*phi/(alpha*dt*bw); %
                    
                    r(n,1)=-cop*(po(n)-pon(n))-cow*(sw(n)-swn(n))+qo(n);
                    r(n+dim,1)=-cwp*(po(n)-pon(n))-cww*(sw(n)-swn(n))+qw(n);
                    
                    jac(n,n)=-cop;
                    jac(n,n+dim)=-cow;
                    
                    jac(n+dim,n)=-cwp;
                    jac(n+dim,n+dim)=-cww;
                    
                    % -x:
                    if (i ~= 1)
                        
                        miuo1=miu_o(po(n-1));
                        miuw1=miu_w(po(n-1));
                        [kro1, kro_s1]=kr_o(sw(n-1));
                        [ krw1, krw_s1]=kr_w(sw(n-1));
                        [bo1, inv_bo_p]=b_o(po(n-1));
                        [bw1, inv_bw_p]=b_w(po(n-1));
                        if po(n)>=po(n-1)
                            to1=beta*0.5*ax*(kx(n)+kx(n-1))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro;
                            tw1=beta*0.5*ax*(kx(n)+kx(n-1))/dx* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw;
                            
                            %kx(n-1)
                            r(n,1)=r(n,1)+to1*(po(n-1)-po(n)-gamma_o*(z(n-1)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n-1)-po(n)-gamma_w*(z(n-1)-z(n)));
                            
                            jac(n, n)=jac(n,n)-to1 ;
                            jac(n,n-1)=to1 ;
                            
                            jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ax*(kx(n)+kx(n-1))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n-1)-po(n)-gamma_o*(z(n-1)-z(n)));
                            % jac(n, n+dim-1)=0;
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1 ;
                            jac(n+dim, n-1)=tw1 ;
                            
                            jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*ax*(kx(n)+kx(n-1))/dx* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n)-po(n-1)-gamma_w*(z(n)-z(n-1)));
                            % jac(n+dim, n+dim-1)=0;
                        else
                            to1=beta*0.5*ax*(kx(n)+kx(n-1))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1;
                            tw1=beta*0.5*ax*(kx(n)+kx(n-1))/dx* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1;
                            
                            r(n,1)=r(n,1)+to1*(po(n-1)-po(n)-gamma_o*(z(n-1)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n-1)-po(n)-gamma_w*(z(n-1)-z(n)));
                            
                            jac(n, n)=jac(n,n)-to1 ;
                            jac(n,n-1)=to1 ;
                            
                            % jac(n, n+dim)=jac(n, n+dim)+0;
                            jac(n, n+dim-1)=beta*0.5*ax*(kx(n)+kx(n-1))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n-1)-po(n)-gamma_o*(z(n-1)-z(n)));
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1 ;
                            jac(n+dim, n-1)=tw1 ;
                            
                            % jac(n+dim, n+dim)=jac(n+dim, n+dim)+0.0;
                            jac(n+dim, n+dim-1)=beta*0.5*ax*(kx(n)+kx(n-1))/dx* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n-1)-po(n)-gamma_w*(z(n-1)-z(n)));
                        end
                        
                    end
                    
                    % +x:
                    if (i ~= nx)
                        miuo1=miu_o(po(n+1));
                        miuw1=miu_w(po(n+1));
                        [kro1, kro_s1]=kr_o(sw(n+1));
                        [krw1, krw_s1]=kr_w(sw(n+1));
                        [bo1, inv_bo_p]=b_o(po(n+1));
                        [bw1, inv_bw_p]=b_w(po(n+1));
                        if po(n)>=po(n+1)
                            to1=beta*0.5*ax*(kx(n)+kx(n+1))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro;
                            tw1=beta*0.5*ax*(kx(n)+kx(n+1))/dx* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw;
                            
                            r(n,1)=r(n,1)+to1*(po(n+1)-po(n)-gamma_o*(z(n+1)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n+1)-po(n)-gamma_w*(z(n+1)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n+1)=to1;
                            
                            jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ax*(kx(n)+kx(n+1))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n+1)-po(n)-gamma_o*(z(n+1)-z(n)));
                            % jac(n, n+dim+1)=0;
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n+1)=tw1;
                            
                            jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*ax*(kx(n)+kx(n+1))/dx* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n+1)-po(n)-gamma_w*(z(n+1)-z(n)));
                            % jac(n+dim, n+dim+1)=0;
                        else
                            to1=beta*0.5*ax*(kx(n)+kx(n+1))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1;
                            tw1=beta*0.5*ax*(kx(n)+kx(n+1))/dx* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1;
                            
                            r(n,1)=r(n,1)+to1*(po(n+1)-po(n)-gamma_o*(z(n+1)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n+1)-po(n)-gamma_w*(z(n+1)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n+1)=to1;
                            
                            % jac(n, n+dim)=jac(n, n+dim)+0;
                            jac(n, n+dim+1)=beta*0.5*ax*(kx(n)+kx(n+1))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n+1)-po(n)-gamma_o*(z(n+1)-z(n)));
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n+1)=tw1;
                            
                            % jac(n+dim, n+dim)=jac(n+dim, n+dim)+0;
                            jac(n+dim, n+dim+1)=beta*0.5*ax*(kx(n)+kx(n+1))/dx* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n+1)-po(n)-gamma_w*(z(n+1)-z(n)));
                        end
                        
                        
                    end
                    
                    % -y:
                    if (j ~= 1)
                        miuo1= miu_o(po(n-nx));
                        miuw1=miu_w(po(n-nx));
                        [kro1, kro_s1]=kr_o(sw(n-nx));
                        [krw1, krw_s1]=kr_w(sw(n-nx));
                        [bo1, inv_bo_p]= b_o(po(n-nx));
                        [bw1, inv_bw_p]= b_w(po(n-nx));
                        if po(n)>=po(n-nx)
                            to1=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro;
                            tw1=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw;
                            
                            r(n,1)=r(n,1)+to1*(po(n-nx)-po(n)-gamma_o*(z(n-nx)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n-nx)-po(n)-gamma_w*(z(n-nx)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n-nx)=to1;
                            
                            jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ay*(ky(n)+ky(n-nx))/dy* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n-nx)-po(n)-gamma_o*(z(n-nx)-z(n)));
                            % jac(n, n+dim-nx)=0;
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n-nx)=tw1;
                            
                            jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*ay*(ky(n)+ky(n-nx))/dy* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n-nx)-po(n)-gamma_w*(z(n-nx)-z(n)));
                            % jac(n+dim, n+dim-nx)=0;
                        else
                            to1=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1;
                            tw1=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1;
                            
                            r(n,1)=r(n,1)+to1*(po(n-nx)-po(n)-gamma_o*(z(n-nx)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n-nx)-po(n)-gamma_w*(z(n-nx)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n-nx)=to1;
                            
                            % jac(n, n+dim)=jac(n, n+dim)+0;
                            jac(n, n+dim-nx)=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n-nx)-po(n)-gamma_o*(z(n-nx)-z(n)));
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n-nx)=tw1;
                            
                            % jac(n+dim, n+dim)=jac(n+dim, n+dim)+0;
                            jac(n+dim, n+dim-nx)=beta*0.5*ay*(ky(n)+ky(n-nx))/dy* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n-nx)-po(n)-gamma_w*(z(n-nx)-z(n)));
                        end
                        
                        
                    end
                    
                    % +y:
                    if (j ~= ny)
                        miuo1= miu_o(po(n+nx));
                        miuw1= miu_w(po(n+nx));
                        [kro1, kro_s1]=kr_o(sw(n+nx));
                        [krw1, krw_s1]=kr_w(sw(n+nx));
                        [bo1, inv_bo_p]=b_o(po(n+nx));
                        [bw1, inv_bw_p]= b_w(po(n+nx));
                        if po(n)>=po(n+nx)
                            to1=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro;
                            tw1=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw;
                            
                            r(n,1)=r(n,1)+to1*(po(n+nx)-po(n)-gamma_o*(z(n+nx)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n+nx)-po(n)-gamma_w*(z(n+nx)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n+nx)=to1;
                            
                            jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ay*(ky(n)+ky(n+nx))/dy* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n+nx)-po(n)-gamma_o*(z(n+nx)-z(n)));
                            % jac(n, n+dim+nx)=0;
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n+nx)=tw1;
                            
                            jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*ay*(ky(n)+ky(n+nx))/dy* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n+nx)-po(n)-gamma_w*(z(n+nx)-z(n)));
                            % jac(n+dim, n+dim+nx)=0;
                        else
                            to1=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1;
                            tw1=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1;
                            
                            r(n,1)=r(n,1)+to1*(po(n+nx)-po(n)-gamma_o*(z(n+nx)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n+nx)-po(n)-gamma_w*(z(n+nx)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n+nx)=to1;
                            
                            % jac(n, n+dim)=jac(n, n+dim)+0;
                            jac(n, n+dim+nx)=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n+nx)-po(n)-gamma_o*(z(n+nx)-z(n)));
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n+nx)=tw1;
                            
                            % jac(n+dim, n+dim)=jac(n+dim, n+dim)+0;
                            jac(n+dim, n+dim+nx)=beta*0.5*ay*(ky(n)+ky(n+nx))/dy* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n+nx)-po(n)-gamma_w*(z(n+nx)-z(n)));
                        end
                        
                    end
                    
                    % -z:
                    if (k ~= 1)
                        miuo1=miu_o(po(n-nx*ny));
                        miuw1=miu_w(po(n-nx*ny));
                        [kro1, kro_s1]=kr_o(sw(n-nx*ny));
                        [krw1, krw_s1]=kr_w(sw(n-nx*ny));
                        [bo1, inv_bo_p]=b_o(po(n-nx*ny));
                        [bw1, inv_bw_p]=b_w(po(n-nx*ny));
                        if po(n)>=po(n-nx*ny)
                            to1=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro;
                            tw1=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw;
                            
                            r(n,1)=r(n,1)+to1*(po(n-nx*ny)-po(n)-gamma_o*(z(n-nx*ny)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n-nx*ny)-po(n)-gamma_w*(z(n-nx*ny)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n-nx*ny)=to1;
                            
                            jac(n, n+dim)=jac(n, n+dim)+beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n-nx*ny)-po(n)-gamma_o*(z(n-nx*ny)-z(n)));
                            % jac(n, n+dim-nx*ny)=0;
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n-nx*ny)=tw1;
                            
                            jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n-nx*ny)-po(n)-gamma_w*(z(n-nx*ny)-z(n)));
                            % jac(n+dim, n+dim-nx*ny)=0;
                        else
                            to1=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))* 0.5*(kro+kro1);
                            tw1=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))* 0.5*(krw+krw1);
                            
                            r(n,1)=r(n,1)+to1*(po(n-nx*ny)-po(n)-gamma_o*(z(n-nx*ny)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n-nx*ny)-po(n)-gamma_w*(z(n-nx*ny)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n-nx*ny)=to1;
                            
                            % jac(n, n+dim)=jac(n, n+dim)+0;
                            jac(n, n+dim-nx*ny)=beta*0.5*az*(kz(n)+kz(n-nx*ny))/dz* ....
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n-nx*ny)-po(n)-gamma_o*(z(n-nx*ny)-z(n)));
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n-nx*ny)=tw1;
                            
                            %jac(n+dim, n+dim)=jac(n+dim, n+dim)+0;
                            jac(n+dim, n+dim-nx*ny)=beta*0.5*az*(kz(n)+kz(n-nx))/dz* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n-nx*ny)-po(n)-gamma_w*(z(n-nx*ny)-z(n)));
                        end
                        
                    end
                    
                    % +z:
                    if (k ~= nz)
                        miuo1=miu_o(po(n+nx*ny));
                        miuw1=miu_w(po(n+nx*ny));
                        [kro1, kro_s1]=kr_o(sw(n+nx*ny));
                        [krw1, krw_s1]=kr_w(sw(n+nx*ny));
                        [bo1, inv_bo_p]=b_o(po(n+nx*ny));
                        [bw1, inv_bw_p]=b_w(po(n+nx*ny));
                        if po(n)>=po(n+nx*ny)
                            to1=beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro;
                            tw1=beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw;
                            
                            r(n,1)=r(n,1)+to1*(po(n+nx*ny)-po(n)-gamma_o*(z(n+nx*ny)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n+nx*ny)-po(n)-gamma_w*(z(n+nx*ny)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n+nx*ny)=to1;
                            
                            jac(n, n+dim)=jac(n, n+dim)+beta*0.5*ax*(kx(n)+kx(n+nx*ny))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s*(po(n+nx*ny)-po(n)-gamma_o*(z(n+nx*ny)-z(n)));
                            % jac(n, n+dim+nx*ny)=0;
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n+nx*ny)=tw1;
                            
                            jac(n+dim, n+dim)=jac(n+dim, n+dim)+beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s*(po(n+nx*ny)-po(n)-gamma_w*(z(n+nx*ny)-z(n)));
                            % jac(n+dim, n+dim+nx*ny)=0;
                        else
                            to1=beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro1;
                            tw1=beta*0.5*az*(kz(n)+kz(n+nx*ny))/dz* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw1;
                            
                            r(n,1)=r(n,1)+to1*(po(n+nx*ny)-po(n)-gamma_o*(z(n+nx*ny)-z(n)));
                            r(n+dim,1)=r(n+dim,1)+tw1*(po(n+nx*ny)-po(n)-gamma_w*(z(n+nx*ny)-z(n)));
                            
                            jac(n, n)=jac(n, n)-to1;
                            jac(n,n+nx*ny)=to1;
                            
                            % jac(n, n+dim)=jac(n, n+dim)+0;
                            jac(n, n+dim+nx*ny)=beta*0.5*ax*(kx(n)+kx(n+nx*ny))/dx* ...
                                0.5*(1/(miuo*bo)+1/(miuo1*bo1))*kro_s1*(po(n+nx*ny)-po(n)-gamma_o*(z(n+nx*ny)-z(n)));
                            
                            jac(n+dim, n)=jac(n+dim, n)-tw1;
                            jac(n+dim, n+nx*ny)=tw1;
                            
                            % jac(n+dim, n+dim)=jac(n+dim, n+dim)+0;
                            jac(n+dim, n+dim+nx*ny)=beta*0.5*az*(kz(n)+kz(n-nx))/dz* ...
                                0.5*(1/(miuw*bw)+1/(miuw1*bw1))*krw_s1*(po(n+nx*ny)-po(n)-gamma_w*(z(n+nx*ny)-z(n)));
                        end
                        
                    end
                    
                end
            end
        end
        
        DX=-inv(jac)*r;

                for p=1:dim
                    po(p)=po(p)+DX(p,1);
                    sw(p)=sw(p)+DX(p+dim,1);
                end
        t
        error= norm(DX)
    end
    pon=po;
    swn=sw;
end
cputime-time