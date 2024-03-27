%function anaSeg()
%%
%clear
thr=0.6;
minSz=15;

seg=tom_mrcread('T27bin2_DNA__DNA_seg_final_cleaned_forangleCalculation.mrc');
seg=seg.Value;

%skel=tom_skeleton3D(seg>thr);


stats = regionprops3(seg>thr,'all');

vol=stats.Volume;
idx=find(vol>minSz);

vtmp=zeros(size(seg));
allIdx=[];
for i=1:length(idx) 
    allIdx=cat(1,allIdx,stats.VoxelIdxList{idx(i)}); 
end
vtmp(allIdx)=1;
% tom_vol2chimera(v);

avgEv=[0 0 0]';
for i=1:length(idx)
   x(i)=stats.Centroid(idx(i),1);
   y(i)=stats.Centroid(idx(i),2);
   z(i)=stats.Centroid(idx(i),3);
   ev=stats.EigenVectors{idx(i)};
   ev1=ev(:,1);
   u(i)=ev1(1);
   v(i)=ev1(2); 
   w(i)=ev1(3); 
   
   if (atan2d(norm(cross([1 0 0],ev1)), dot([1 0 0],ev1))>90)
        avgEv=avgEv+(ev1.*-1);
   else
       avgEv=avgEv+ev1;
   end
end

avgEv=avgEv./length(idx);


figure; hold on; axis image;
quiver3(x,y,z,v,u,w,'color','#1A75BC'); 
%%

% quiver3(Clustrx,Clustry,Clustrz,Clustrv,Clustru,Clustrw,'color','#DB80AD');
%%
title('principle axis of filaments');

for i=1:length(idx)
     ev=stats.EigenVectors{idx(i)};
     ev1=ev(:,1);
     angle(i) = atan2d(norm(cross(avgEv,ev1)), dot(avgEv,ev1));
     if (angle(i)>90)
        angle(i)=180-angle(i);
     end
end
[N X] = hist(angle,100);
f = fit(X.',N.','gauss1')
figure; histogram(angle,100);
hold on; plot(f)
% [N2 X2] = hist(Clustrangle,100);
% f2 = fit(X2.',N2.','gauss1')
% histogram(Clustrangle,100);
% hold on; plot(f2)
title('histogram of angles between average principle axis and individual principle axis');







