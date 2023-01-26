function [x,y,z] = L2xyz(L,delta)
% covert the size of a box to the xyz coordinates in 1D

Lx = L(1); Ly = L(2); Lz = L(3);
dx = delta(1); dy = delta(2); dz = delta(3);

Nx = Lx/dx; Ny = Ly/dy; Nz = Lz/dz;
x           = dx* (-Nx/2+0.5:1:Nx/2-0.5);      % 1D axis in x
y           = dy* (-Ny/2+0.5:1:Ny/2-0.5);      % 1D axis in y
z           = dz* (-Nz/2+0.5:1:Nz/2-0.5);      % 1D axis in y

end

