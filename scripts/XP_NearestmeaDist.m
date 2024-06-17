function meaDist()
 
starName='run_data.star';
nrBin=400; % bins for hist
distCut=0;%  %in voxel
pixsize = 0.352;% in nm
filtStarName='xxxxx.star';% star file from Relion
pathToTOM=''; %abs path to tomFolde from nemotoc project

%% Code
addpath(genpath(pathToTOM));
st=tom_starread(starName);
dSt=tom_extractData(st);

for i=1:length(st) 
    motl(8,i)=st(i).rlnCoordinateX; 
    motl(9,i)=st(i).rlnCoordinateY;   
    motl(10,i)=st(i).rlnCoordinateZ; 
    motl(5,i)=dSt.label.tomoID(i); 
end

idRemove=[];
labelu=unique(dSt.label.tomoID);
for i=1:length(labelu) 
    idx=find(dSt.label.tomoID==labelu(i));  
    for ii=1:length(idx)
        d=sg_pairwise_dist(motl(8:10,idx(ii)),motl(8:10,idx)); 
        d(ii)=100000;
        allD(idx(ii))=min(d);
        if (min(d)<distCut)
            indNeigh=[find(d<distCut)];
            idRemove=cat(1,idRemove,idx(indNeigh));
        end
    end 
end
allD = allD*pixsize;
[N X] = hist(allD,nrBin);
f = fit(X.',N.','gauss1')
figure; hist(allD,nrBin);
hold on; plot(f)
disp([num2str(length(idRemove)) ' particles got removed by dist']);
idGood=ismember(1:length(st),idRemove)==0;
stF=st(idGood);
stF(1).Header=st(1).Header;
disp(['Writing cleaned star to: ' filtStarName]);
tom_starwrite(filtStarName,stF);



function D = sg_pairwise_dist(X, Y)
% A function to determine the pairwise distance between two n-dimensional
% arrays. I found it after googling "matlab pairwise distance". It should
% be substantially faster than the matlab function pdist.
%
% A fix had to be implemented as the bsxfun changed in the matlab 8.x
% and produces rounding errors close to zero. This results in some numbers 
% becoming negative, making the distance array a complex array. This fix
% rounds all negative values to zero. 
%
% WW 02-2016


% Squares of the distances
sq_dist = (bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y));

% Find negative values
neg = sq_dist < 0;
% Set negative values to zero
sq_dist(neg) = 0;

D = sq_dist.^0.5;

