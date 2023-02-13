function [RI] = MakeSphereInRandMed(rad, n, L, delta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a single (stratified) sphere imbeded in a random medium
%
% Lang Wang, 2023
% -------------------------------------------------------------------------
%                              INPUTS
%
% rad            -> radii of the stritified sphere, from center to edge
% n              -> the RI of the sphere and the background
%
% -------------------------------------------------------------------------
%                              OUTPUTS
%
% dCsdO          -> the normalized RCS
% ang            -> angles of the detectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_imm = n(end);

[x,y,z] = L2xyz(L,delta);
[X, Y, Z]=meshgrid(x,y,z);

Nx = size(x,2); Ny = size(y,2); Nz = size(z,2);
RI=n_imm*ones(Nx,Ny,Nz);
for ii = 1:size(rad,2)
dielectricSphere=X.^2+Y.^2+Z.^2<=rad(ii)^2;
RI=RI + double(dielectricSphere)*(n(ii)-n(ii+1));
end

figure
subplot(1,2,1)
imagesc(x,y,RI(:,:,end/2));
xlabel('x(\lambda)')
ylabel('y(\lambda)')
axis square
colorbar
set(gca,'YDir','normal')
title('xy plane')
subplot(1,2,2)
imagesc(x,z,squeeze(RI(:,end/2,:)));
xlabel('x(\lambda)')
ylabel('z(\lambda)')
axis square
colorbar
set(gca,'YDir','normal')
title('xz plane')

end

