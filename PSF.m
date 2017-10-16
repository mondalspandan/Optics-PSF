a=127;   %z-sensitivity
b=127;   %r-sensitivity

alpha=45*pi/180;
delta=5*pi/180;			%10*pi/180 (used 10 for the plots in the report)
lambda=976*10^-9;
k=2*pi/lambda;

zmax=5*10^-6;
rmax=zmax*sin(alpha);

z=linspace(-zmax,zmax,a);
r=linspace(-rmax,rmax,b);

[Z,R]=meshgrid(z,r);

u=k.*z.*(sin(alpha))^2;
v=k.*r.*sin(alpha);

[U,V]=meshgrid(u,v);

phi=0*pi/180;
phicirc=45*pi/180;

A=39.83;

I0AE1 = zeros(b,a);
I1AE1 = zeros(b,a);
I2AE1 = zeros(b,a);
I0AE2 = zeros(b,a);
I1AE2 = zeros(b,a);
I2AE2 = zeros(b,a);
I0BW = zeros(b,a);
I1BW = zeros(b,a);
I2BW = zeros(b,a);

for j = 1:b
    for i = 1:a
        fnc0 = @(th) (cos(th)).^.5.*sin(th).*(1+cos(th)).*besselj(0,v(j).*sin(th)./sin(alpha)).*exp(1i*u(i).*cos(th)./(sin(alpha))^2);
        I0AE1(j,i)= quad(fnc0,0,delta);
        I0AE2(j,i)= quad(fnc0,alpha-delta,alpha);
        I0BW(j,i)= quad(fnc0,0,alpha);
        
        fnc1 = @(th) (cos(th)).^.5.*(sin(th)).^2.*besselj(1,v(j).*sin(th)./sin(alpha)).*exp(1i*u(i).*cos(th)./(sin(alpha))^2);
        I1AE1(j,i)= quad(fnc1,0,delta);
        I1AE2(j,i)= quad(fnc1,alpha-delta,alpha);
        I1BW(j,i)= quad(fnc1,0,alpha);
        
        fnc2 = @(th) (cos(th)).^.5.*sin(th).*(1-cos(th)).*besselj(2,v(j).*sin(th)./sin(alpha)).*exp(1i*u(i).*cos(th)./(sin(alpha))^2);
        I2AE1(j,i)= quad(fnc2,0,delta);
        I2AE1(j,i)= quad(fnc2,alpha-delta,alpha);
        I2BW(j,i)= quad(fnc2,0,alpha);
    end
end

hcpBW=(A^2/(8*pi))^2.* (abs(I0BW).^2+2.*abs(I1BW).^2+abs(I2BW).^2).^3;

hlpBW=(A^2/(8*pi))^2.* (  abs(I0BW).^2  +  4.*abs(I1BW).^2.*cos(phi).^2  +  abs(I2BW).^2  +  2.*cos(2.*phi).*real(I0BW.*conj(I2BW))  ).^2    .*   (abs(I0BW).^2  +  2.*abs(I1BW).^2  +  abs(I2BW).^2);

ex = -1i*A.*(I0AE1+I2AE1.*cos(2.*phi))-1i*A.*(I0AE2+I2AE2.*cos(2.*phi));
ey = -1i*A.*I2AE1.*sin(2.*phi)-1i*A.*I2AE2.*sin(2.*phi);
ez = -2*A.*I1AE1.*cos(phi)-2*A.*I1AE2.*cos(phi);    
E=(abs(ex).^2+abs(ey).^2+abs(ez).^2).^.5;

exBW = -1i*A.*(I0BW+I2BW.*cos(2.*phi));
eyBW = -1i*A.*I2BW.*sin(2.*phi);
ezBW = -2*A.*I1BW*cos(phi);    
EBW=(abs(exBW).^2+abs(eyBW).^2+abs(ezBW).^2).^.5;

excirc=-1i*A.*(I0AE1+I2AE1.*cos(2*phicirc))-1i*A.*(I0AE2+I2AE2.*cos(2*phicirc));
eycirc=-1i*A.*I2AE1.*sin(2*phicirc)-1i*A.*I2AE2.*sin(2*phicirc);
ezcirc=-2*A.*I1AE1.*cos(phicirc)-2*A.*I1AE2.*cos(phicirc);
Ecirc=(abs(excirc).^2+abs(eycirc).^2+abs(ezcirc).^2).^.5;

excircBW=-1i*A.*(I0BW+I2BW.*cos(2*phicirc));
eycircBW=-1i*A.*I2BW.*sin(2*phicirc);
ezcircBW=-2*A.*I1BW.*cos(phicirc);
EcircBW=(abs(excircBW).^2+abs(eycircBW).^2+abs(ezcircBW).^2).^.5;

hlpAE = (abs(ex).^2+abs(ey).^2+abs(ez).^2).^2  +  (abs(I0AE1+I0AE2).^2+2.*abs(I1AE1+I1AE2).^2+abs(I2AE1+I2AE2).^2);

hcpAE=(abs(excirc).^2+abs(eycirc).^2+abs(ezcirc).^2).^2  +  (abs(I0AE1+I0AE2).^2+2.*abs(I1AE1+I1AE2).^2+abs(I2AE1+I2AE2).^2);

figure
set(gcf, 'Position', [750, 750, 750, 750])
subplot(2,2,1)
plot3(U,V,hcpAE)
title('AE CP')
subplot(2,2,2)
plot3(U,V,hlpAE)
title('AE LP')
subplot(2,2,3)
plot3(U,V,hcpBW)
title('BW CP')
subplot(2,2,4)
plot3(U,V,hcpBW)
title('BW LP')

figure
set(gcf, 'Position', [750, 750, 750, 750])
subplot(2,2,1)
contourf(U,V,hcpAE)
title('AE CP')
subplot(2,2,2)
contourf(U,V,hlpAE)
title('AE LP')
subplot(2,2,3)
contourf(U,V,hcpBW)
title('BW CP')
subplot(2,2,4)
contourf(U,V,hlpBW)
title('BW LP')

figure
set(gcf, 'Position', [750, 750, 750, 750])
subplot(2,2,1)
plot(V,hcpAE(:,(a+1)/2))
title('AE CP')
subplot(2,2,2)
plot(V,hlpAE(:,(a+1)/2))
title('AE LP')
subplot(2,2,3)
plot(V,hcpBW(:,(a+1)/2))
title('BW CP')
subplot(2,2,4)
plot(V,hcpBW(:,(a+1)/2))
title('BW LP')

figure
set(gcf, 'Position', [750, 750, 750, 750])
subplot(2,2,1)
contourf(U,V,Ecirc)
title('AE CP')
subplot(2,2,2)
contourf(U,V,E)
title('AE LP')
subplot(2,2,3)
contourf(U,V,EcircBW)
title('BW CP')
subplot(2,2,4)
contourf(U,V,EBW)
title('BW LP')

rcpAE=hcpAE(:,(a+1)/2);
HM=rcpAE((b+1)/2)/2;
pts=[];
for rval=1:b-1
    if (rcpAE(rval)-HM)*(rcpAE(rval+1)-HM)<=0
        pts=[pts rval];
    end
end
FWHMforAEcp=(pts(2)-pts(1))*2*rmax/b

rcpBW=hcpBW(:,(a+1)/2);
HM=rcpBW((b+1)/2)/2;
pts=[];
for rval=1:b-1
    if (rcpBW(rval)-HM)*(rcpBW(rval+1)-HM)<=0
        pts=[pts rval];
    end
end
FWHMforBWcp=(pts(2)-pts(1))*2*rmax/b

rlpAE=hlpAE(:,(a+1)/2);
HM=rlpAE((b+1)/2)/2;
pts=[];
for rval=1:b-1
    if (rlpAE(rval)-HM)*(rlpAE(rval+1)-HM)<=0
        pts=[pts rval];
    end
end
FWHMforAElp=(pts(2)-pts(1))*2*rmax/b

rlpBW=hlpBW(:,(a+1)/2);
HM=rlpBW((b+1)/2)/2;
pts=[];
for rval=1:b-1
    if (rlpBW(rval)-HM)*(rlpBW(rval+1)-HM)<=0
        pts=[pts rval];
    end
end
FWHMforBWlp=(pts(2)-pts(1))*2*rmax/b
