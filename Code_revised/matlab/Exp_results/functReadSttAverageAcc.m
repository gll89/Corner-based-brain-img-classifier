function [acc_arr, pre_arr, rec_arr, fscore_arr] = functReadSttAverageAcc(path, f, folder, thSubstr)
% path is /9result/Classifierx/resultAll
suf_txt = ['\*Fold' num2str(f) '.txt'];
txtNames = dir([path, suf_txt]);
% txtNames = dir([path, '\*Fold3.txt']);
% txtNames = {txtNames(:).name};
leng = length(txtNames);
acc_arr = zeros(10,10);
pre_arr = zeros(10,10);  %[];
rec_arr = zeros(10,10);  %[];
fscore_arr = zeros(10,10); %[];

kArr = 3:2:21;
for i = 1:leng   %read all the statistics text
    txtName = txtNames(i).name;
    txtPath = fullfile(path, txtName);
    C = importdata(txtPath);
    dataC = C{1, 1};
    strArr = strsplit(dataC);
    lengD = [3:3:21];  %the indices of numbers in C
    dataD = strArr(lengD);  % the values of acc, rec, pre, spe, npv,mcc...
    for j = 1:length(dataD)
        tmp = dataD(j);
        tmp = cell2mat(tmp);
        dataTmp = str2num(tmp);
        data(1, j) = dataTmp;
    end
    acc = data(1,1);
    pre = data(1,2);    
    rec = data(1,3);
    fscore = 2*pre*rec/(pre+rec); 
  
    addNum = 5;
    lengtxtName = length(txtName);
    if strcmp(folder , 'CT')
        if lengtxtName == 19+addNum
            th = txtName(10);
            k = txtName(13);
        else if lengtxtName == 20+addNum
                if ~strcmp(txtName(11), 'T') 
                    th = txtName(10:11);
                    k = txtName(14);
                else 
                    th = txtName(10);
                    k = txtName(13:14);
                end
            else if lengtxtName == 21+addNum
            th = txtName(10:11);
            k = txtName(14:15);
                end
            end
        end
        th = str2num(th);
        k = str2num(k);
        k = find(kArr == k);  
        
        acc_arr(k, th)= acc;
        pre_arr(k, th) = pre;
        rec_arr(k, th) = rec;
        fscore_arr(k, th) = fscore;
%         resultAll(k, th-thSubstr) = dataArr;
    end
    
    if strcmp(folder, 'MRI')
        %%----HarrisMRI
        if lengtxtName == 19+addNum
            th = txtName(10);
            k = txtName(13);
        else if lengtxtName == 21+addNum
            th = txtName(10:11);
            k = txtName(14:15);
            else if lengtxtName == 20+addNum
                    if ~strcmp(txtName(11), 'T')
                    th = txtName(10:11);
                    k = txtName(14);
                else 
                    th = txtName(10);
                    k = txtName(13:14);
                end
                end
            end
        end
        th = str2num(th);
        k = str2num(k);
        k = find(kArr == k);  
        th = th - thSubstr;
        
        acc_arr(k, th)= acc;
        pre_arr(k, th) = pre;
        rec_arr(k, th) = rec;
        fscore_arr(k, th) = fscore;

    end
    
end
end


