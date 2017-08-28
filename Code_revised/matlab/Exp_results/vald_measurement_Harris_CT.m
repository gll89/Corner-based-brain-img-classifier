clear all, clc, close all;

%% -----------MRI-------------------
folder =  'CT';  
th = 1:10; 
thSubstr = 0;  % used for adjustment of th value between the real value and the  value
th_range =  0.001:0.001:0.01;
thmin = 1; thmax = 10;
kmin=3; kmax=21;
k = 3:2:21;
path = 'E:\GLL-BMC\BMC\Image_revised\HarrisCT\9result\measurement';
savefile='acc_arr_harris_ct.mat';
% acc_arr_3d = [];
% for f = 1:5
%     [acc_arr, pre_arr, rec_arr, fscore_arr] = functReadSttAverageAcc(path, f, folder, thSubstr);   %row is k, and column is th
%     maxV = max(acc_arr(:));
%     [maxIndX, maxIndY] = find(acc_arr == maxV);
%     acc_arr_3d(:, :, f) = acc_arr;
% end
% save(savefile, 'acc_arr_3d');
load(savefile, 'acc_arr_3d');
figure(3),
%% ======subplot======
for f = 1:5
    subplot(3, 2,f);
    
    acc_mat = acc_arr_3d(:,:,f);
    surf(k, th, acc_mat, 'LineWidth', 0.3);
    axis tight
    colormap summer ;  %spring;
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

