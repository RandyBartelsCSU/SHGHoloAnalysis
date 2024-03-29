clear variables; close all; clc; addpath(genpath('../Functions'));

ps          = .105*2;                 % pixel size (x,y,z) in object space (microns)
lambda      = 0.5*2;                  % central wavelength (microns)
R_test = 20*lambda;
NA          = 0;                  % numerical aperture of imaging and detection lens
n_imm       = 1;                % refractive index of immersion media
nsphere=1.2;
rad = 2*lambda;
n=[nsphere,n_imm];
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = 2^9;                  % lateral pixel dimension 
L = [N*ps, N*ps, N*ps];
delta = [ps, ps, ps];

[x,y,z] = L2xyz(L,delta);
[X,Z]=meshgrid(x,z);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);
dGk = 1;
%z_inc = z(1)-ps;
Eps=0.5/lambda^2;

if NA==0
U_inp=ones(N,N);
end

ord = 1;

E_MSR=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,ord,Eps,dGk,'Vol');
E_MLR=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MLB=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017

dR = ps;
[D_R_MSR,TH] = RCS(E_MSR,X,Z,k,R_test,dR);
[D_R_MLR,TH] = RCS(E_MLR,X,Z,k,R_test,dR);
[D_R_MLB,TH] = RCS(E_MLB,X,Z,k,R_test,dR);

%[D_R_MSR,TH] = RCS_line(E_MSR,x,z,k,R_test,dR);
%[D_R_MLR,TH] = RCS_line(E_MLR,x,z,k,R_test,dR);
%[D_R_MLB,TH] = RCS_line(E_MLB,x,z,k,R_test,dR);

[dCsdO, ang]= MieRCS(rad, n, lambda);

c=299792458;
[an,bn,RCSTheta,ETheta] = mieHKURCS(rad,c/lambda,n_imm^2,1,nsphere^2,1,40,TH);

figure
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
pbaspect([3 1 1])
