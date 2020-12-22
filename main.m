close all;
clear all;
clc;

%-------
% Paths.
%-------
addpath('modelspecific');
addpath('multigs');

%-------------------
% Compile Mex files.
%-------------------
cd multigs;
if exist('computeIntersection','file')~=3
    % 编译c文件
    mex computeIntersection.c;
end
cd ..;

%----------------------
% Setup VLFeat toolbox.
%----------------------
% cd vlfeat-0.9.14/toolbox;
% feval('vl_setup');
% cd ../..;
run('./vlfeat-0.9.14/toolbox/vl_setup');

%---------------------------------------------
% Check if we are already running in parallel.
%---------------------------------------------
% poolsize = matlabpool('size');
% if poolsize == 0 %if not, we attempt to do it:
%     matlabpool open;
% end
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers
end

%-------------------------
% User defined parameters.
%-------------------------
% Global model specific function handlers.
clear global;
global fitfn resfn degenfn psize numpar
fitfn = 'homography_fit';
resfn = 'homography_res';
degenfn = 'homography_degen';
psize   = 4;
numpar  = 9;

M     = 500;  % Number of hypotheses for RANSAC.
thr   = 0.1;  % RANSAC threshold.

scale = 1;    % Scale of input images (maybe for large images you would like to use a smaller scale).

%------------------
% Images to stitch.
%------------------
example_dir = '/home/lynx/study/stitch/dataset/test/'
path1 = [example_dir '1.jpg'];
path2 = [example_dir '2.jpg'];

%-------------
% Read images.
%-------------
fprintf('Read images and SIFT matching\n');tic;
fprintf('> Reading images...');tic;
img1 = imresize(imread(sprintf('%s',path1)),scale);
img2 = imresize(imread(sprintf('%s',path2)),scale);
fprintf('done (%fs)\n',toc);

%--------------------------------------
% SIFT keypoint detection and matching.
%--------------------------------------
fprintf('  Keypoint detection and matching...');tic;
[ kp1,ds1 ] = vl_sift(single(rgb2gray(img1)),'PeakThresh', 0,'edgethresh',500);
[ kp2,ds2 ] = vl_sift(single(rgb2gray(img2)),'PeakThresh', 0,'edgethresh',500);
matches   = vl_ubcmatch(ds1,ds2);
fprintf('done (%fs)\n',toc);

% Normalise point distribution.
% 归一化
fprintf('  Normalising point distribution...');tic;
data_orig = [ kp1(1:2,matches(1,:)) ; ones(1,size(matches,2)) ; kp2(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];
[ dat_norm_img1,T1 ] = normalise2dpts(data_orig(1:3,:));
[ dat_norm_img2,T2 ] = normalise2dpts(data_orig(4:6,:));
data_norm = [ dat_norm_img1 ; dat_norm_img2 ];
fprintf('done (%fs)\n',toc);

%-----------------
% Outlier removal.
%-----------------
fprintf('Outlier removal\n');tic;
% Multi-GS
rng(0);
[ ~,res,~,~ ] = multigsSampling(100,data_norm,M,10);
con = sum(res<=thr);
[ ~, maxinx ] = max(con);
inliers = find(res(:,maxinx)<=thr);

% if size(img1,1) == size(img2,1)
%     % Show results of RANSAC.
%     fprintf('  Showing results of RANSAC...');tic;
%     figure;
%     imshow([img1 img2]);
%     hold on;
%     plot(data_orig(1,:),data_orig(2,:),'ro','LineWidth',2);
%     plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'ro','LineWidth',2);
%     for i=1:length(inliers)
%         plot(data_orig(1,inliers(i)),data_orig(2,inliers(i)),'go','LineWidth',2);
%         plot(data_orig(4,inliers(i))+size(img1,2),data_orig(5,inliers(i)),'go','LineWidth',2);
%         plot([data_orig(1,inliers(i)) data_orig(4,inliers(i))+size(img1,2)],[data_orig(2,inliers(i)) data_orig(5,inliers(i))],'g-');
%     end
%     title('Ransac''s results');
%     fprintf('done (%fs)\n',toc);
% end
