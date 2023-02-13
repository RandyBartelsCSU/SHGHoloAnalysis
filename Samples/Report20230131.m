clear variables; close all; clc; addpath(genpath('../Functions'));

ps_array = [0.105, 0.07, 0.08, 0.09,0.12];

for ii = 1%:size(ps_array,2)
    ii

ps          = ps_array(ii);%.105;                 % pixel size (x,y,z) in object space (microns)
lambda      = 0.5;                  % central wavelength (microns)
NA          = 0;                  % numerical aperture of imaging and detection lens
n_imm       = 1;                % refractive index of immersion media
nsphere=1.02;
n=[nsphere,n_imm];
k0=(2*pi)/lambda;
k=k0*n_imm;
%N_round = round(2^9*.105/ps);
%if rem(N_round,2)==1
%    N_round=N_round+1;
%end
N           = [2^9, 2^9, 2^9];%[2^9, 2^9, 2^9];     [N_round, N_round, N_round]             % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];
dGk = 1; % 0 to 1

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

rad = 2*lambda;

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);

%z_inc = z(1)-ps;
Eps=0.5/lambda^2;

if NA==0
U_inp=ones(N(1),N(2));
end

ord = 1;

% Gk clip

bound = dGk * ( max(max(sqrt(fxx.^2 + fyy.^2))) - n_imm/lambda );
SquareRt=@(a) abs(real(sqrt(a)))+1i*abs(imag(sqrt(a))); % making sure imaginary and real part is pos otherwise will get gain or backward propagating field
prop_phs= 1i*2*pi*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2));%+1i*eps %Add small absorption term to avoid indeterminates in the angular greens function
prop=@(z) exp(prop_phs*z);
%Mask=(n_imm/lambda)^2>1.01*((fxx.^2+fyy.^2));
AG= ((-1i.*exp(prop_phs.*ps)./(4.*pi.*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2)+1i*Eps)))); % Angular Greens function
AG(fxx.^2 + fyy.^2 > (n_imm/lambda + bound)^2) = 0;
AG = fftshift(AG);
AG_slide = squeeze(AG(:,end/2))';

plot(fftshift(fx),abs(AG_slide));

E_MSR=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,ord,Eps,dGk,'Vol');
E_MLR=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MLB=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017

% Last xy plane

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
[~,~,RCSTheta1D,ETheta1D] = mieHKURCS(rad,c/lambda,n_imm^2,1,nsphere^2,1,40,theta1D);
ETheta = reshape(ETheta1D, [size(theta,1) size(theta,2)]);

R = sqrt(X.^2 + Y.^2 + z(end)^2);
E_plane = ETheta./( exp(1i*k*z(end))/z(end) ).*( exp(1i*k*R)./R );
E_plane = E_plane./E_plane(end/2,end/2);

% RCS

[X,Z]=meshgrid(x,z);

R_test = 20*lambda;
dR = ps;
[D_R_MSR,~] = RCS(E_MSR,X,Z,k,R_test,dR);
[D_R_MLR,~] = RCS(E_MLR,X,Z,k,R_test,dR);
[D_R_MLB,TH] = RCS(E_MLB,X,Z,k,R_test,dR);

[dCsdO, ang] = MieRCS(rad, n, lambda);

c=299792458;
[an,bn,RCSTheta,ETheta] = mieHKURCS(rad,c/lambda,n_imm^2,1,nsphere^2,1,40,TH);

% xz plane

E_MSR_xz = squeeze(E_MSR(end/2,:,:));
E_MLR_xz = squeeze(E_MLR(end/2,:,:));
E_MLB_xz = squeeze(E_MLB(end/2,:,:));

% Plot
L_plot = 10;

figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'papertype','A4');
f_title = sprintf('n = %1.2f, r = %1.2f \\lambda, \\Delta G_k = %1.2f, [ps] = %1.2f \\lambda, \\epsilon = %1.2f \\lambda^{-2}', nsphere, rad./lambda, dGk, ps/lambda, Eps*lambda^2);
sgtitle(f_title)

% object
subplot(4,7,1)
imagesc(x./lambda,y./lambda,RI(:,:,end/2));
xlabel('x(\lambda)')
ylabel('y(\lambda)')
axis square
colorbar
set(gca,'YDir','normal')
title('xy plane')
xlim([-L_plot/2 L_plot/2])
ylim([-L_plot/2 L_plot/2])
subplot(4,7,2)
imagesc(x./lambda,z./lambda,squeeze(RI(:,end/2,:)));
xlabel('x(\lambda)')
ylabel('z(\lambda)')
axis square
colorbar
set(gca,'YDir','normal')
title('xz plane')
xlim([-L_plot/2 L_plot/2])
ylim([-L_plot/2 L_plot/2])

% object
subplot(4,7,3)
plot(fftshift(fx),abs(AG_slide));

% xz field
cmax = max(max(abs(E_MSR_xz)));
subplot(4,7,15)
imagesc(x./lambda,y./lambda,abs(E_MSR_xz));
xlabel('x(\lambda)')
ylabel('z(\lambda)')
clim([0 cmax]);
axis square
colorbar
set(gca,'YDir','normal')
title('MSR |E_{tot}|')
subplot(4,7,22)
imagesc(x./lambda,y./lambda,angle(E_MSR_xz));
xlabel('x(\lambda)')
ylabel('z(\lambda)')
axis square
colorbar
set(gca,'YDir','normal')
title('MSR \angle E_{tot}')

subplot(4,7,16)
imagesc(x./lambda,y./lambda,abs(E_MLR_xz));
xlabel('x(\lambda)')
ylabel('z(\lambda)')
clim([0 cmax]);
axis square
colorbar
set(gca,'YDir','normal')
title('MLR |E_{tot}|')
subplot(4,7,23)
imagesc(x./lambda,y./lambda,angle(E_MLR_xz));
xlabel('x(\lambda)')
ylabel('z(\lambda)')
axis square
colorbar
set(gca,'YDir','normal')
title('MLR \angle E_{tot}')

subplot(4,7,17)
imagesc(x./lambda,y./lambda,abs(E_MLB_xz));
xlabel('x(\lambda)')
ylabel('z(\lambda)')
clim([0 cmax]);
axis square
colorbar
set(gca,'YDir','normal')
title('MLB |E_{tot}|')
subplot(4,7,24)
imagesc(x./lambda,y./lambda,angle(E_MLB_xz));
xlabel('x(\lambda)')
ylabel('z(\lambda)')
axis square
colorbar
set(gca,'YDir','normal')
title('MLB \angle E_{tot}')

% RCS
cmax = max(max(abs(E_sca_MSR)));
subplot(4,7,[8 9 10])
hold on
plot(180./pi*TH,10*log10(D_R_MSR),'LineWidth',4);
plot(180./pi*TH,10*log10(D_R_MLR),'LineStyle',"--",'LineWidth',4);
plot(180./pi*TH,10*log10(D_R_MLB),'LineStyle',":",'LineWidth',4);
plot(ang, 10*log10(dCsdO),'LineStyle',"-.",'LineWidth',4);
plot(180./pi*TH,10*log10(RCSTheta),'LineStyle',"-",'LineWidth',4);
xlim([0 80])
ylim([-40 0])
xlabel('Scattering angle (deg)')
ylabel('Normalized RCS (dB)')
legend('MSR','MLR','MLB','Mie MatScat','Mie HKU')
%pbaspect([3 1 1])

subplot(4,7,7)
imagesc(x./lambda,y./lambda,abs(E_plane))
clim([0 cmax]);
axis square
colorbar
title('ground truth |Mie|')

subplot(4,7,14)
imagesc(x./lambda,y./lambda,angle(E_plane))
axis square
colorbar
title('ground truth \angle Mie')

subplot(4,7,4)
imagesc(x./lambda,y./lambda,abs(E_sca_MSR))
axis square
clim([0 cmax]);
colorbar
title('MSR |E_{sca}|')

subplot(4,7,11)
imagesc(x./lambda,y./lambda,angle(E_sca_MSR))
axis square
colorbar
title('MSR \angle E_{sca}')

subplot(4,7,18)
imagesc(x./lambda,y./lambda,abs(E_sca_MSR-E_plane))
clim([0 cmax]);
axis square
colorbar
title('MSR |E_{sca}-E_{Mie}|')

subplot(4,7,25)
imagesc(x./lambda,y./lambda,angle(E_sca_MSR-E_plane))
axis square
colorbar
title('\angle E_{sca}-E_{Mie}')

subplot(4,7,5)
imagesc(x./lambda,y./lambda,abs(E_sca_MLR))
axis square
clim([0 cmax]);
colorbar
title('MLR |E_{sca}|')

subplot(4,7,12)
imagesc(x./lambda,y./lambda,angle(E_sca_MLR))
axis square
colorbar
title('\angle E_{sca}')

subplot(4,7,19)
imagesc(x./lambda,y./lambda,abs(E_sca_MLR-E_plane))
clim([0 cmax]);
axis square
colorbar
title('MLR |E_{sca}-E_{Mie}|')

subplot(4,7,26)
imagesc(x./lambda,y./lambda,angle(E_sca_MLR-E_plane))
axis square
colorbar
title('\angle E_{sca}-E_{Mie}')

subplot(4,7,6)
imagesc(x./lambda,y./lambda,abs(E_sca_MLB))
axis square
clim([0 cmax]);
colorbar
title('MLB |E_{sca}|')

subplot(4,7,13)
imagesc(x./lambda,y./lambda,angle(E_sca_MLB))
axis square
colorbar
title('\angle E_{sca}')

subplot(4,7,20)
imagesc(x./lambda,y./lambda,abs(E_sca_MLB-E_plane))
clim([0 cmax]);
axis square
colorbar
title('MLB |E_{sca}-E_{Mie}|')

subplot(4,7,27)
imagesc(x./lambda,y./lambda,angle(E_sca_MLB-E_plane))
axis square
colorbar
title('\angle E_{sca}-E_{Mie}')

saveas(gcf,[pwd sprintf('/Figures/ps%d.jpg',ii)]);

end