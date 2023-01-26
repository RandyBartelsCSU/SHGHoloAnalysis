function [R fieldoutsize]=GenerateReflectionMatrix(Centx,Centy,SidelobeFilter,HoloFilter,domainswitch,Size,filename)



holoinfo=h5info(filename,'/Epi/Hologram');
count=holoinfo.ChunkSize;
Epiinfo=h5info(filename,'/Epi');
datainfo=Epiinfo.Datasets.Dataspace;
datasize=datainfo.Size;
rootinfo=h5info(filename);
rootinfo.Attributes.Value   % display data comments

ps=(1/40)/83.5; %mm
ps=ps*1000; %um             % pixel size from callibration image
Fs=1/ps;

Nx=datasize(2);
Ny=datasize(1);  
x=ps*[-Nx/2:Nx/2-1];
y=ps*[-Ny/2:Ny/2-1];
dFx=Fs/Nx;
dfxs         = dFx;%1/(N*ps);             % Fourier spacing 
fxs          = dfxs*[-Nx/2:Nx/2-1];         % 1D axis in fx
dFy=Fs/Ny;
dfys         = dFy;%1/(N*ps);             % Fourier spacing 
fys          = dfys*[-Ny/2:Ny/2-1];         % 1D axis in fy
[Xs Ys]=meshgrid(x,y);
[fxxs fyys]=meshgrid(fxs,fys);

SidelobeFilter=abs(ifft2(ifftshift(fftshift(fft2(SidelobeFilter)).*exp(1i*2*pi*((Centx)*dfxs*Xs+Centy*dfys*Ys))))); % Shift sidelobe filter to location of centroid


switch Size
    case 'Full'

R=zeros(datasize(1)*datasize(2),datasize(3));
switch domainswitch
    case 'field'
upd = textprogressbar(datasize(3), 'barlength', 20, ...
                         'updatestep', 1, ...
                         'startmsg', 'Generating reflection matrix... ',...
                         'endmsg', ' Finally!', ...
                         'showbar', true, ...
                         'showremtime', true, ...
                         'showactualnum', true, ...
                         'barsymbol', '+', ...
                         'emptybarsymbol', '-');
for ii=1:datasize(3)
    start=[1 1 ii];

    ACHolo=h5read(filename,'/Epi/Hologram',start,count)-h5read(filename,'/Epi/Reference',start,count)-h5read(filename,'/Epi/Signal',start,count);
    ACHolo=ACHolo.*(HoloFilter);
    fftholo=(fftshift(fft2((ACHolo)))).*SidelobeFilter; % apply filter
    % Appy phase ramp to shift to center
    field=ifft2(ifftshift(fftholo)).*exp(1i*2*pi*((Centx)*dfxs*Xs+Centy*dfys*Ys)); %

    R(:,ii)=field(:);
    upd(ii);

end
fieldoutsize=size(field);
    case 'kspace'
upd = textprogressbar(datasize(3), 'barlength', 20, ...
                         'updatestep', 1, ...
                         'startmsg', 'Generating reflection matrix... ',...
                         'endmsg', ' Finally!', ...
                         'showbar', true, ...
                         'showremtime', true, ...
                         'showactualnum', true, ...
                         'barsymbol', '+', ...
                         'emptybarsymbol', '-');
for ii=1:datasize(3)
    start=[1 1 ii];

    ACHolo=h5read(filename,'/Epi/Hologram',start,count)-h5read(filename,'/Epi/Reference',start,count)-h5read(filename,'/Epi/Signal',start,count);
    ACHolo=ACHolo.*(HoloFilter);
    fftholo=(fftshift(fft2((ACHolo)))).*SidelobeFilter; % apply filter
    % Appy phase ramp to shift to center
    field=ifft2(ifftshift(fftholo)).*exp(1i*2*pi*((Centx)*dfxs*Xs+Centy*dfys*Ys)); %
    fieldk=fftshift(fft2(field));

    R(:,ii)=fieldk(:);
    upd(ii);

end
fieldoutsize=size(fieldk);
otherwise
        disp('Choose between field and kspace for the option domainswitch')
end

    case 'Trimmed'


        
switch domainswitch
    case 'field'
upd = textprogressbar(datasize(3), 'barlength', 20, ...
                         'updatestep', 1, ...
                         'startmsg', 'Generating reflection matrix... ',...
                         'endmsg', ' Finally!', ...
                         'showbar', true, ...
                         'showremtime', true, ...
                         'showactualnum', true, ...
                         'barsymbol', '+', ...
                         'emptybarsymbol', '-');

ff=figure; 
        start=[1 1 round((datasize(3)-1)/2+1)];

    Holo=h5read(filename,'/Epi/Hologram',start,count);
imagesc(Holo)
colormap Turbo
%caxis([0 1e4])
disp('Choose FOV...')
[xfove yfove]=ginput(2);
yfove=round(yfove);
xfove=round(xfove);
close(ff);

ff=figure; 
imagesc(Holo(yfove(1):yfove(2),xfove(1):xfove(2))) 
colormap Turbo
pause(1)
close(ff);

% set the FOV 
xMine = xfove(1);
xMaxe = xfove(2);
yMine = yfove(1);
yMaxe = yfove(2);






R=zeros((yMaxe-yMine+1)*(xMaxe-xMine+1),datasize(3));
fprintf('\n Output field size is %4.2f by %8.3f \n',yMaxe-yMine+1,yMaxe-yMine+1)
%disp(['The output size for the fields is ' num2str(yMaxe-yMine+1) 'by' num2str(xMaxe-xMine+1)])
for ii=1:datasize(3)
    start=[1 1 ii];

    ACHolo=h5read(filename,'/Epi/Hologram',start,count)-h5read(filename,'/Epi/Reference',start,count)-h5read(filename,'/Epi/Signal',start,count);
    ACHolo=ACHolo.*(HoloFilter);
    fftholo=(fftshift(fft2((ACHolo)))).*SidelobeFilter; % apply filter
    % Appy phase ramp to shift to center
    field=ifft2(ifftshift(fftholo)).*exp(1i*2*pi*((Centx)*dfxs*Xs+Centy*dfys*Ys)); %
    fieldt=field(yfove(1):yfove(2),xfove(1):xfove(2));

    R(:,ii)=fieldt(:);
    upd(ii);

end
fprintf('\n Output field size is %4.2f by %8.3f \n',yMaxe-yMine+1,yMaxe-yMine+1)
fieldoutsize=size(fieldt);
    case 'kspace'
fybnd=max(find(fys<-.69 & fys>-.71));  % These bounds set by the NA and wavelength
fxbnd=max(find(fxs<-.69 & fxs>-.71));
    
R=zeros(length(fys(fybnd:end-fybnd))*length(fxs(fxbnd:end-fxbnd)),datasize(3));
fprintf('\n Output field size is %4.2f by %8.3f \n',length(fys(fybnd:end-fybnd)),length(fxs(fxbnd:end-fxbnd)))

upd = textprogressbar(datasize(3), 'barlength', 20, ...
                         'updatestep', 1, ...
                         'startmsg', 'Generating reflection matrix... ',...
                         'endmsg', ' Finally!', ...
                         'showbar', true, ...
                         'showremtime', true, ...
                         'showactualnum', true, ...
                         'barsymbol', '+', ...
                         'emptybarsymbol', '-');

for ii=1:datasize(3)
    start=[1 1 ii];

    ACHolo=h5read(filename,'/Epi/Hologram',start,count)-h5read(filename,'/Epi/Reference',start,count)-h5read(filename,'/Epi/Signal',start,count);
    ACHolo=ACHolo.*(HoloFilter);
    fftholo=(fftshift(fft2((ACHolo)))).*SidelobeFilter; % apply filter
    % Appy phase ramp to shift to center
    field=ifft2(ifftshift(fftholo)).*exp(1i*2*pi*((Centx)*dfxs*Xs+Centy*dfys*Ys)); %
    fieldk=fftshift(fft2(field));
    
    fieldkt=fieldk(fybnd:end-fybnd,fxbnd:end-fxbnd);

    R(:,ii)=fieldkt(:);
    upd(ii);

end

fprintf('\n Output field size is %4.2f by %8.3f \n',length(fys(fybnd:end-fybnd)),length(fxs(fxbnd:end-fxbnd)))
fieldoutsize=size(fieldkt);
otherwise
        disp('Choose between field and kspace for the option domainswitch')
end

end

end