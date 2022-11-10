function [FilterHolo FilterSidelobe]=FilterGen(filename,ps,NA,lambda)


% This function generates filters both for the hologram in R-space and for
% the filtering of the sidelobe containing field information in k-space for
% the purpososes of off axis holography. The filter in R-space is designed
% to mitigate ringing artifacts when taking fourier transforms of the
% hologram due to non-zero values near the edge of the field of view. The
% filter in the fourier domain is used subsequently to select and filter
% one of the sidelobes to reconstruct the given field

% Inputs:

% ps:       Calibrated pixel size corresponding to a camera image in
%           R-space. Units: microns

% NA:       Numerical aperture of collection objective

% Lambda:   Wavelength of collected light





holoinfo=h5info(filename,'/Epi/Hologram');
count=holoinfo.ChunkSize;
Epiinfo=h5info(filename,'/Epi');
datainfo=Epiinfo.Datasets.Dataspace;
datasize=datainfo.Size;

Nx=datasize(2);
Ny=datasize(1);


x=linspace(-1,1,Nx);
y=linspace(-1,1,Ny);

[X Y]=meshgrid(x,y);

sig2=.7;
FilterHolo=exp(-((X).^2/sig2).^4).*exp(-((Y).^2/sig2).^4);
figure; 


imagesc(FilterHolo)
daspect([1 1 1])
axis off
title('Filter for R-space Hologram')
figure;
tiledlayout(1,2)
nexttile
plot(FilterHolo(round(end/2),:))
title('y-dir lineout')
nexttile
plot(FilterHolo(:,round(end/2)))
title('x-dir lineout')
% pause(3)
% close(ff)

% ps=(1/40)/83.5; %mm
% ps=ps*1000;     %um


Fs=1/ps;
dFx=Fs/Nx;
dfxs= dFx;%1/(N*ps);             % Fourier spacing 
fxs= dfxs*[-Nx/2:Nx/2-1];% 1D axis in fx
dFy=Fs/Ny;
dfys= dFy;%1/(N*ps);             % Fourier spacing 
fys= dfys*[-Ny/2:Ny/2-1];% 1D axis in fy

[fxxs fyys]=meshgrid(fxs,fys);

% NA=.16;
% lambda=1.030;

sig=(NA/lambda);
%Filter=exp(-(((X-s(2)/2+cent(1)).^2+(Y-s(1)/2).^2)/sig).^4);
FilterSidelobe=exp(-(((fxxs).^2+(fyys).^2)/sig).^4);

figure; 

imagesc(fxs,fys,FilterSidelobe)
title('Sidelobe filter')
xlabel('f_x')
ylabel('f_y')
daspect([1 1 1])
hold on
circle(0,0,sig)
legend('NA cutoff')

% pause(3)
% close(ff)


end

function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'k--');
end