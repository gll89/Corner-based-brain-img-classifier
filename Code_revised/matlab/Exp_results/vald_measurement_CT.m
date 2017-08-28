clear all, clc, close all;

%% --------CTImgs folder
folder =  'CT';  
th = 10:19; 
thSubstr = 9;  % used for adjustment of th value between the real value and the  value
th_range = 0.01:0.01:0.1;
thmin = 10; thmax = 19;
kmin=3; kmax=21;

%% 
savefile='acc_arr_ct.mat';
k = 3:2:21;
path = 'E:\GLL-BMC\BMC\Image_revised\CTImgs\9result\meansurement_all_folds';
acc_arr_3d = [];
% for f = 1:5
%     [acc_arr, pre_arr, rec_arr, fscore_arr] = functReadSttAverageAcc(path, f, folder, thSubstr);   %row is k, and column is th
%     maxV = max(acc_arr(:));
%     [maxIndX, maxIndY] = find(acc_arr == maxV);
%     acc_arr_3d(:, :, f) = acc_arr;
% end
% save(savefile, 'acc_arr_3d');
load(savefile, 'acc_arr_3d');
%% ======subplot======
for f = 1:5
    subplot(3, 2,f);
    % figure(2),
    acc_mat = acc_arr_3d(:,:,f);
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

    maxV_tmp = max(acc_mat(:));
    [xind, yind] = find( acc_mat == maxV_tmp);
    maxV = zeros(size(xind));
    maxV(:) = maxV_tmp;
    h = scatter3(k(yind), th(xind), maxV, 'filled');
    h.SizeData = 120;
    hold off;
    disp(['f=', num2str(f), ', maxV=', num2str(maxV(1)), ', th=', num2str(xind(1)), ', k=', num2str(yind(1))]);
end

% %% the average acc
% acc_mat = mean(acc_arr_3d, 3);
% subplot(3, 2,6);
% surf(k, th, acc_mat, 'LineWidth', 0.3);
% axis tight
% colormap cool;
% colorbar
% axis([kmin, kmax, thmin, thmax,  min(acc_mat(:)), max(acc_mat(:))]);
% xlabel('K');
% ylabel('\theta');
% zlabel('Accuracy');
% set(gca, 'XTick', [3 5 7 9 11 13 15 17 19 21]);
% set(gca, 'XTicklabel', [3 5 7 9 11 13 15 17 19 21]);
% set(gca, 'YTick', th );
% set(gca, 'YTicklabel', th_range);
% hold on,
% 
% maxV_tmp = max(acc_mat(:));
% [xind, yind] = find( acc_mat == maxV_tmp);
% maxV = zeros(size(xind));
% maxV(:) = maxV_tmp;
% h = scatter3(k(yind), th(xind), maxV, 'filled');
% h.SizeData = 120;
% hold off;


