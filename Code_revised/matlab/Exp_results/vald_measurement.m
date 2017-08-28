clear all, clc, close all;
root_path = 'E:\Download\Dropbox\BMC Bioinformatics-2017\First-round-review\Image_revised';
extension_path = '9result\meansurement';

%% produce the 3D figure of validation data under different thresholds of th and K

%% --------CT  Our proposed method---------
folder =  'CT';  
th = 10:19; 
thSubstr = 9;  % used for adjustment of th value between the real value and the  value
th_range = 0.01:0.01:0.1;
thmin = 10; thmax = 19;
kmin=3; kmax=21;

% read experimental results from text files
savefile='acc_arr_ct.mat';
% load(savefile, 'acc_arr_3d');
k = 3:2:21;
path = fullfile(root_path, 'CTImgs', extension_path);
f=3;
[acc_mat, pre_arr, rec_arr, fscore_arr] = functReadSttAverageAcc(path, f, folder, thSubstr);   %row is k, and column is th
save(savefile, 'acc_mat');
% ======subplot figures======
figure(1),
subplot(2, 1, 1);
surf(k, th, acc_mat, 'LineWidth', 0.3);
axis tight
colormap cool;
colorbar
axis([kmin, kmax, thmin, thmax,  min(acc_mat(:)), max(acc_mat(:))]);
xlabel('K');
ylabel('\theta');
zlabel('Accuracy');
set(gca, 'XTick', [3 5 7 9 11 13 15 17 19 21]);
set(gca, 'XTicklabel', [3 5 7 9 11 13 15 17 19 21]);
set(gca, 'YTick', th );
set(gca, 'YTicklabel', th_range);
hold on,
mt(1) = title('(a)');

maxV_tmp = max(acc_mat(:));
[xind, yind] = find( acc_mat == maxV_tmp);
maxV = zeros(size(xind));
maxV(:) = maxV_tmp;
% h = scatter3(k(yind), th(xind), maxV, 'filled', 'markerfacecolor',[1 0 0]);
% h.SizeData = 120;
% hold off;

%% -----------MRI  Our proposed method-------------------
folder =  'MRI';  
th = 3:12; 
thSubstr = 2;  % used for adjustment of th value between the real value and the  value
th_range = 0.03:0.01:0.12;
thmin = 3; thmax = 12;
kmin=3; kmax=21;
k = 3:2:21;
path = fullfile(root_path, 'MRIImgs', extension_path);
acc_arr_3d = [];
f=3;
savefile='acc_arr_mri.mat';
% load(savefile, 'acc_arr_3d');
[acc_mat, pre_arr, rec_arr, fscore_arr] = functReadSttAverageAcc(path, f, folder, thSubstr);   %row is k, and column is th
save(savefile, 'acc_arr_3d');
% ======subplot======
subplot(2, 1, 2);
surf(k, th, acc_mat, 'LineWidth', 0.3);
axis tight
%     colormap winter;
colormap cool;
colorbar
axis([kmin, kmax, thmin, thmax,  min(acc_mat(:)), max(acc_mat(:))]);
xlabel('K');
ylabel('\theta');
zlabel('Accuracy');
set(gca, 'XTick', [3 5 7 9 11 13 15 17 19 21]);
set(gca, 'XTicklabel', [3 5 7 9 11 13 15 17 19 21]);
set(gca, 'YTick', th );
set(gca, 'YTicklabel', th_range);
hold on,

maxV_tmp = max(acc_mat(:));
[xind, yind] = find( acc_mat == maxV_tmp);
xind = xind(2);
yind = yind(2);
maxV = zeros(size(xind));
maxV(:) = maxV_tmp;
% h = scatter3(k(yind), th(xind), maxV,  'marker','x', 'markerfacecolor',[1 0 0]);
% h.SizeData = 120;
% hold off;
% mt(2) = 
title('(b)');
text(5, 0.4, 'Bottom title')


%% --------Harris CT---------
folder =  'CT';  
th = 10:19; 
thSubstr = 9;  % used for adjustment of th value between the real value and the  value
th_range = 0.001:0.001:0.01;
thmin = 10; thmax = 19;
kmin=3; kmax=21;

% read experimental results from text files
savefile='acc_arr_ct.mat';
% load(savefile, 'acc_arr_3d');
k = 3:2:21;
path = fullfile(root_path, 'HarrisCT', extension_path);
f=3;
[acc_mat, pre_arr, rec_arr, fscore_arr] = functReadSttAverageAcc(path, f, folder, thSubstr);   %row is k, and column is th
save(savefile, 'acc_mat');
% ======subplot figures======
figure(2),
subplot(2, 1, 1);
surf(k, th, acc_mat, 'LineWidth', 0.3);
axis tight
colormap default; 
colorbar
axis([kmin, kmax, thmin, thmax,  min(acc_mat(:)), max(acc_mat(:))]);
xlabel('K');
ylabel('\theta');
zlabel('Accuracy');
set(gca, 'XTick', [3 5 7 9 11 13 15 17 19 21]);
set(gca, 'XTicklabel', [3 5 7 9 11 13 15 17 19 21]);
set(gca, 'YTick', th );
set(gca, 'YTicklabel', th_range);
set(gca, 'fontsize', 8);
hold on,
title('(a)');

% maxV_tmp = max(acc_mat(:));
% [xind, yind] = find( acc_mat == maxV_tmp);
% maxV = zeros(size(xind));
% maxV(:) = maxV_tmp;
% h = scatter3(k(yind), th(xind), maxV, 'filled', 'markerfacecolor',[1 0 0]);
% h.SizeData = 120;
% hold off;

subplot(2, 1, 2);
surf(k, th, acc_mat, 'LineWidth', 0.3);
axis tight
colormap default;
colorbar
axis([kmin, kmax, thmin, thmax,  min(acc_mat(:)), max(acc_mat(:))]);
xlabel('K');
ylabel('\theta');
zlabel('Accuracy');
set(gca, 'XTick', [3 5 7 9 11 13 15 17 19 21]);
set(gca, 'XTicklabel', [3 5 7 9 11 13 15 17 19 21]);
set(gca, 'YTick', th );
set(gca, 'YTicklabel', th_range);
set(gca, 'fontsize', 8);
hold on,
title('(b)');
