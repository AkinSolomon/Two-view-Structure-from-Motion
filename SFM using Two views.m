clc
clear all
%Akinlawon Solomon
%run('C:\Users\Akinlawon\Documents\MATLAB\vlfeat-0.9.20\toolbox\vl_setup')
%need vl_sift package to run code

%% Initialization Step

imfile = dir('images2/');
imgs = {imfile(~[imfile.isdir]).name};
randimg = [407:6:462];

for i = 1:size(randimg,2)
    imgss{1,i} = imread(imgs{1,randimg(i)});
end

peakthresh = 1;
matchthresh = 2;
I1col=imgss{1,1};
I1=single(rgb2gray(I1col));
[f1, d1] = vl_sift(I1, 'PeakThresh', 1);
pics=10;
match_points=zeros(2,1500,pics-1);
imgfeat{1,1}=f1(1:2,:);
imgdesc{1,1}=d1;


for j=1:length(randimg)-1
    I2 =imgss{1,j+1};
    I2=single(rgb2gray(I2));
    
    [f2, d2] = vl_sift(I2, 'PeakThresh', peakthresh);
    [matches, scores] = vl_ubcmatch(d1, d2, matchthresh);
    num_matches = length(matches);
    imgfeat{1,j+1} = f2(1:2,:);
    imgdesc{1,j+1} = d2;
    match_points(1:2,1:num_matches,j)=matches;
    matchesPerImage(j) = num_matches;
end

perm = randperm(size(f1,2));
sel = perm(1:50);
figure(1)
imshow(I1)
hold on
h1 = vl_plotframe(f1(:,sel));
h2 = vl_plotframe(f1(:,sel));
set(h1,'color','k','linewidth',3);
set(h2,'color','y','linewidth',2);
title('SIFT frames at a subset of points of interest')

figure(2)
imshow(I1)
h3 = vl_plotsiftdescriptor(d1(:,sel),f1(:,sel));
set(h3,'color','g');
title('SIFT descriptors at a subset of the points of interest')
%% Find intersection of features
indxset = intersect(match_points(1,:,1), match_points(1,:,2));

for k=3:pics-1
    indxset = intersect(indxset, match_points(1,:,k));
end

indxset=indxset(2:end);
intersectmatch(1,1:length(indxset)) = indxset;

for k=2:pics
    for m=1:length(indxset)
        ind = find(match_points(1,:,k-1)==indxset(m));
        intersectmatch(k,m) = match_points(2,ind,k-1);
    end
end

p=zeros(3*pics,length(indxset));
one=ones(1,length(indxset));


for m=1:pics
    f=imgfeat{1,m};
    p(3*(m-1)+1:3*(m-1)+3,:)=[f(1:2,intersectmatch(m,:));one];
end

plotmatch(imgss{1,1},imgss{1,i},f1,f2,indxset);

[~,ind]=unique(p(1,:));
ind=sort(ind);
p=p(:,ind);
Wnorm=[];


for m=1:pics
    temp = p(3*(m-1)+1:3*(m-1)+3,:);
    tmp = normalise2dpts(temp);
    Wnorm=[Wnorm; tmp];
end

W=p;
W2=Wnorm;

%
%% Calibration matrix
x1 = W2(1:3,:);
cx =(max(x1(1,:))- min(x1(1,:)))/2;
cy =(max(x1(2,:))- min(x1(2,:)))/2;
fx = 2*cx-min(x1(2,:));
fy = fx;
K = [fx 0 cx;0 fy cy;0 0 1];


clear P M
[P M] = triangulate(W2,K,size(W2,2),1,10);
P=P(:,P(3,:)<50);

figure;
scatter3(NaN,NaN,NaN);
hold on;

for k=1:length(P)
    
   W_point = W(1:2,k);
   orX = round(W_point(1));
   orY = round(W_point(2));
   orRGB = I1col(orY, orX,:);
   scatter3(P(1,k), P(2,k), P(3,k), 'MarkerEdgeColor', double(squeeze(orRGB)')/255, 'MarkerFaceColor', double(squeeze(orRGB)')/255);
   
end

hold off
title('SFM using two other views')
axis equal

%%
