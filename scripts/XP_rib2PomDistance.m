% calculate distances of ribosomes to Pom cluster
% Xu peng @MPI 2023_01_12
% plot ribosome density map back
% Pom cluster segmentation from EMAN

%segFolder = '/fs/gpfs03/lv03/pool/pool-plitzko/Peng_Xu/FromMJ/tomography/t8/seg/amira';
%cd(segFolder)

INMrec = 'T27bin4_cluster__cluster_seg.mrc';% segmentation of pom cluster
chromrec = 'plotback_clean_2nd.mrc';% plot back of ribosomes

whichINM = INMrec; % from which INM do you want to measure the distance?


ps=0.352; % pixelsize in nm
bin=4; % binning factor
resample=1; % rate that you want to resample your data at, 1 is fine
gaussians=5; % how many Gaussian distributions do you want to fit to the data?
histBins = 1000; % number of bins for the histogram



% measure distances from INM
% create a distance matrix 
mat = tom_mrcread(whichINM);mat = mat.Value;
%tom_volxyz (mat);

% binarize the segmentation, rarely necessary
%mat=double(mat);
%mat(find(mat<0.9))=0;
%mat(find(mat>0.9))=1;

% generate 'distance matrix'
mat1 = bwdist(mat);
%tom_volxyz (mat1);


%% 
chrom=tom_mrcread(chromrec);chrom = chrom.Value;
%tom_volxyz (chrom);

% binarize if necessary
%chrom=double(chrom);
%chrom(find(chrom<0.9))=0;
%chrom(find(chrom>0.9))=1;

% transfrom to 1D vector containnig only 1s and resample data
v=find(chrom>0);
out=randperm(size(v,1));
vr=v(out);
vr_resample=vr(1:resample:end);
coor=zeros(size(vr_resample,1),3);
%[x y z] = ind2sub ([size(chrom,1) size(chrom,2) size(chrom,3)],vr_resample);
[coor(:,1) coor(:,2) coor(:,3)] = ind2sub ([size(chrom,1) size(chrom,2) size(chrom,3)],vr_resample);

% get distances in pixels
for i =1:size(coor,1)
     coor(i,4)=mat1(round(coor(i,1)),round(coor(i,2)),round(coor(i,3)));
end
 
% get distances in nm
coor(:,5)=coor(:,4)*ps*bin;


% fit gauss and get mean and std
[N X]=hist(coor(:,5),histBins);

f = fit(X.',N.',['gauss' num2str(gaussians)])

figure; hist(coor(:,5),histBins);
hold on; plot(f)

