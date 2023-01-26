function [] = DORT(R,fieldoutsize)
Nx=fieldoutsize(1);
Ny=fieldoutsize(2);
pix = 1;
[x, y, X, Y] = Make_Mesh(Nx,Ny,pix);

N_O_rec = 24;
[G_rec,Lambda,H_rec] = svds(R,N_O_rec);
%G_rec_sum = zeros(N_OBJ, N_OBJ);
%H_rec_sum = zeros(Nx, Ny);
%GH_rec_sum = zeros(N_OBJ, N_OBJ);

for ii=1:N_O_rec
    ii
    %G_rec_col = reshape(G_rec(:,ii),[N_CCD, N_CCD]);
    H_rec_row = reshape(H_rec(:,ii)',[Ny, Nx]);
    figure(3)
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(4,6,ii)
    imagesc(x, y, abs(H_rec_row));
    set(gca,'YDir','normal')
    colormap hot
    colorbar
    daspect([1 1 1])
    caxis([0 0.1])
    titles = sprintf('H(%d)',ii);
    title(titles)
end
end

