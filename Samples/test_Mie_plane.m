clear variables; close all; clc; addpath(genpath('../Functions'));

ps          = .105;                 % pixel size (x,y,z) in object space (microns)
lambda      = 0.5;                  % central wavelength (microns)
NA          = 0;                  % numerical aperture of imaging and detection lens
n_imm       = 1;                % refractive index of immersion media
nsphere=1.2;
n=[nsphere,n_imm];
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = [2^9, 2^9, 2^9];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

dGk = 1;
rad = 2*lambda;

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);

%z_inc = z(1)-ps;
Eps=1;

if NA==0
U_inp=ones(N(1),N(2));
end

ord = 1;

E_MSR=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,ord,Eps,dGk,'Vol');
E_MLR=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MLB=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017

U_inp_end = exp(1i*k*(z(end)-z(1)))*U_inp;

[E_tot_MSR,E_sca_MSR] = tot2sca(E_MSR,U_inp_end);
[E_tot_MLR,E_sca_MLR] = tot2sca(E_MLR,U_inp_end);
[E_tot_MLB,E_sca_MLB] = tot2sca(E_MLB,U_inp_end);

figure
cmax = max(max(abs(E_tot_MSR)));
subplot(2,3,1)
imagesc(x,y,abs(E_tot_MSR))
axis square
clim([0 cmax]);
colorbar
title('MSR |E_{tot}|')

subplot(2,3,4)
imagesc(x,y,abs(E_sca_MSR))
clim([0 cmax]);
axis square
colorbar
title('MSR |E_{sca}|')

subplot(2,3,2)
imagesc(x,y,abs(E_tot_MLR))
axis square
clim([0 cmax]);
colorbar
title('MLR |E_{tot}|')

subplot(2,3,5)
imagesc(x,y,abs(E_sca_MLR))
clim([0 cmax]);
axis square
colorbar
title('MLR |E_{sca}|')

subplot(2,3,3)
imagesc(x,y,abs(E_tot_MLB))
axis square
clim([0 cmax]);
colorbar
title('MLB |E_{tot}|')

subplot(2,3,6)
imagesc(x,y,abs(E_sca_MLB))
clim([0 cmax]);
axis square
colorbar
title('MLB |E_{sca}|')

E_sca_MSR = E_sca_MSR./E_sca_MSR(end/2,end/2);
E_sca_MLR = E_sca_MLR./E_sca_MLR(end/2,end/2);
E_sca_MLB = E_sca_MLB./E_sca_MLB(end/2,end/2);

theta = atan( sqrt(X.^2 + Y.^2)./z(end) );
theta1D = reshape(theta, [1 size(theta,1)*size(theta,2)]);

c=299792458;
[an,bn,RCSTheta1D,ETheta1D] = mieHKURCS(rad,c/lambda,n_imm^2,1,nsphere^2,1,40,theta1D);
ETheta = reshape(ETheta1D, [size(theta,1) size(theta,2)]);

R = sqrt(X.^2 + Y.^2 + z(end)^2);
E_plane = ETheta./( exp(1i*k*z(end))/z(end) ).*( exp(1i*k*R)./R );
E_plane = E_plane./E_plane(end/2,end/2);

figure
cmax = max(max(abs(E_sca_MSR)));

subplot(4,4,4)
imagesc(x,y,abs(E_plane))
clim([0 cmax]);
axis square
colorbar
title('ground truth |Mie|')

subplot(4,4,8)
imagesc(x,y,angle(E_plane))
axis square
colorbar
title('ground truth \angle Mie')

subplot(4,4,1)
imagesc(x,y,abs(E_sca_MSR))
axis square
clim([0 cmax]);
colorbar
title('MSR |E_{sca}|')

subplot(4,4,5)
imagesc(x,y,angle(E_sca_MSR))
axis square
colorbar
title('MSR \angle E_{sca}')

subplot(4,4,9)
imagesc(x,y,abs(E_sca_MSR-E_plane))
clim([0 cmax]);
axis square
colorbar
title('MSR |E_{sca}-E_{Mie}|')

subplot(4,4,13)
imagesc(x,y,angle(E_sca_MSR-E_plane))
axis square
colorbar
title('\angle E_{sca}-E_{Mie}')

subplot(4,4,2)
imagesc(x,y,abs(E_sca_MLR))
axis square
clim([0 cmax]);
colorbar
title('MLR |E_{sca}|')

subplot(4,4,6)
imagesc(x,y,angle(E_sca_MLR))
axis square
colorbar
title('\angle E_{sca}')

subplot(4,4,10)
imagesc(x,y,abs(E_sca_MLR-E_plane))
clim([0 cmax]);
axis square
colorbar
title('MLR |E_{sca}-E_{Mie}|')

subplot(4,4,14)
imagesc(x,y,angle(E_sca_MLR-E_plane))
axis square
colorbar
title('\angle E_{sca}-E_{Mie}')

subplot(4,4,3)
imagesc(x,y,abs(E_sca_MLB))
axis square
clim([0 cmax]);
colorbar
title('MLB |E_{sca}|')

subplot(4,4,7)
imagesc(x,y,angle(E_sca_MLB))
axis square
colorbar
title('\angle E_{sca}')

subplot(4,4,11)
imagesc(x,y,abs(E_sca_MLB-E_plane))
clim([0 cmax]);
axis square
colorbar
title('MLB |E_{sca}-E_{Mie}|')

subplot(4,4,15)
imagesc(x,y,angle(E_sca_MLB-E_plane))
axis square
colorbar
title('\angle E_{sca}-E_{Mie}')

