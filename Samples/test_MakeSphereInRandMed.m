clear variables; close all; clc; addpath(genpath('../Functions'));

rad = [1,2,3];
n = [4, 3, 2, 1];
L = [7, 7, 7];
delta = [0.1, 0.1, 0.1];
[x,y,z] = L2xyz(L,delta);
k0 = 1;
RI = MakeSphereInRandMed(rad, n, L, delta);

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
