clear all, clc, close all;

%% --------CTImgs folder
folder =  'MRI';  
th = 3:12; 
thSubstr = 2;  % used for adjustment of th value between the real value and the  value
th_range = 0.03:0.01:0.12;
thmin = 3; thmax = 12;
kmin=3; kmax=21;

%% 

k = 3:2:21;
path = 'E:\GLL-BMC\BMC\Image_revised\MRIImgs\9result\measurement_all_6';
acc_arr_3d = [];
for f = 1:5
    [acc_arr, pre_arr, rec_arr, fscore_arr] = functReadSttAverageAcc(path, f, thSubstr);   %row is k, and column is th
    maxV = max(acc_arr(:));
    [maxIndX, maxIndY] = find(acc_arr == maxV);
    acc_arr_3d(:, :, f) = acc_arr;
end
% savefile='acc_arr_mri_1.mat';
% load(savefile, 'acc_arr_3d_1');
%% ======subplot======
figure(3),
for f = 1:5
    subplot(3, 2,f);
    % figure(2),
    acc_mat = acc_arr_3d(:,:,f);
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
    maxV = zeros(size(xind));
    maxV(:) = maxV_tmp;
    h = scatter3(k(yind), th(xind), maxV, 'filled');
    h.SizeData = 120;
    hold off;
end
% the average acc
acc_mat = mean(acc_arr_3d, 3);
subplot(3, 2,6);
surf(k, th, acc_mat, 'LineWidth', 0.3);
axis tight
% colormap winter;
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
maxV = zeros(size(xind));
maxV(:) = maxV_tmp;
h = scatter3(k(yind), th(xind), maxV, 'filled');
h.SizeData = 120;
hold off;


