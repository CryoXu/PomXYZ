%% calculte filaments coordinates to center nearest distance
% XuPeng @ MPI Martinsreid 2024_01_10

%% find the center coordinates
clear
% 
segmentation = tom_mrcread(' ')
segmentation = segmentation.Value;
segmentation(find(segmentation>0.03))=1;% bindary cutoff here = 0.03
segmentation(find(segmentation<=0.03))=0;
%tom_mrcwrite(segmentation,'name','test0.02.mrc');
pixsize = 0.293*4;% in nm
%%
% 
[nonZeroY, nonZeroX, nonZeroZ] = ind2sub(size(segmentation), find(segmentation));

% 
minX = min(nonZeroX);
maxX = max(nonZeroX);

% 
minY = min(nonZeroY);
maxY = max(nonZeroY);

% 
minZ = min(nonZeroZ);
maxZ = max(nonZeroZ);



%% calculate segmentaion to center distance


% 
targetX = round((minX+maxX)/2);
targetY = round((minY+maxY)/2);
targetZ = round((minZ+maxZ)/2);

% 
[nonZeroY, nonZeroX, nonZeroZ] = ind2sub(size(segmentation), find(segmentation));

% 
distances = pixsize*sqrt((nonZeroX - targetX).^2 + (nonZeroY - targetY).^2 + (nonZeroZ - targetZ).^2);

% 
figure;
histogram(distances, 'BinWidth', 1); % 
xlabel('Distance to Point ');
ylabel('Frequency');
title('Distance Distribution to Point');

