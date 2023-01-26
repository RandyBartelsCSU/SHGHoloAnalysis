function [fx,fy] = L2fxfy(L,delta)
% convert L to Fourier 2D spatical frequency fx, fy

%[x,y,z] = L2xyz(L,delta);
deltaf = 1./L;
Lf = 1./delta;
[fx,fy,fz] = L2xyz(Lf,deltaf);

fx         = ifftshift(fx);       % FFT shifting Fourier axes
fy         = ifftshift(fy);       % FFT shifting Fourier ax

end

