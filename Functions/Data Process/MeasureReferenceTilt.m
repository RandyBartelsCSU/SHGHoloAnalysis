function [Centx Centy]=MeasureReferenceTilt(filename)

holoinfo=h5info(filename,'/Epi/Hologram');
count=holoinfo.ChunkSize;
Epiinfo=h5info(filename,'/Epi');
datainfo=Epiinfo.Datasets.Dataspace;
datasize=datainfo.Size;

if mod(datasize(3),2)==0
    warning('The scan path should have an odd number of points so that it will pass through the point [0 0].')
    DClocation=round((datasize(3)-1)/2+1);
else

DClocation=(datasize(3)-1)/2+1;
end


start=[1 1 DClocation];

ACHolo=h5read(filename,'/Epi/Hologram',start,count)-h5read(filename,'/Epi/Reference',start,count)-h5read(filename,'/Epi/Signal',start,count);





fftholo=abs(fftshift(fft2(ACHolo)));
ff=figure; 
imagesc(10*log10(fftholo))
title('Select two points which define a box around sideband...')
[xn yn]=ginput(2);
xn=round(xn);
yn=round(yn);
area=pi*(min([(xn(2)-xn(1)) (yn(2)-yn(1))])/2)^2 % Estimate of blob area
area=area*.5;
thresholdValue=mean(10*log10((fftholo(yn(1):yn(2),xn(1):xn(2)))),"all");

close(ff);

binaryImage = 10*log10(fftholo) > thresholdValue;
binaryImage = imfill(binaryImage, 4,'holes');


%subplot(3, 3, 2);
ff=figure; 
tiledlayout(1,2)
nexttile
histogram(10*log10(fftholo),500)
title('Histogram of ACHolo in k-space');
%xlim([0 grayLevels(end)]); % Scale x axis manually.
grid on;

hold on;
maxYValue = ylim;
line([thresholdValue, thresholdValue], maxYValue, 'Color', 'r');
nexttile
imagesc(binaryImage)
pause(1)
close(ff);

fftholo=(fftshift(fft2((ACHolo))));
    %thresholdValue=mean(10*log10((fftholo(yn(1):yn(2),xn(1):xn(2)))),"all");
    binaryImage = 10*log10(abs(fftholo)) > thresholdValue;
    binaryImage = imfill(binaryImage, 4,'holes');

[labeledImage, numberOfBlobs] = bwlabel(binaryImage, 4); 


props = regionprops(labeledImage, 10*log10(abs(fftholo)),'Area','BoundingBox','Centroid');

allBlobAreas = [props.Area];
area=pi*(min([(xn(2)-xn(1)) (yn(2)-yn(1))])/2)^2; % Estimate of blob area
area=area*.5;
allowableAreaIndexes = allBlobAreas > area; % Take the small objects.

keeperIndexes = find(allowableAreaIndexes);
[N, M] = size(fftholo);
Centy=(N+1)/2-props(keeperIndexes(1)).Centroid(2);%-74.1376;
Centx=(M+1)/2-props(keeperIndexes(1)).Centroid(1);%466.6813;


figure; imagesc(10*log10(abs(fftholo)))
daspect([1 1 1])
hold on
plot(props(keeperIndexes(1)).Centroid(1),props(keeperIndexes(1)).Centroid(2),'k*','LineWidth',3)
axis off
colormap hot
legend('Centroid Location')


end


