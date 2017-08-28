root_path = 'E:\GLL-BMC\BMC\Image_revised';

%% =============CT============
method_folder = 'CTImgs';
funct_folder = '9resultTest';
path = fullfile(root_path, method_folder, funct_folder, 'time');
filename = 'ValdRuntime0.02Th3NNFold1.txt';
filepath = fullfile(path, filename);
data = importdata(filepath);
data = data{1};
% class(data)
strArr = strsplit(data);
inds = 4:5:length(strArr);
times = strArr(inds);
img_num = length(times);
time_cell = times(img_num);
time_char = time_cell{1};
time = time_char(1:length(time_char)-2);
time_ave = 1.0 * time*60/img_num

%% =============HarrisCT============
method_folder = 'HarrisCT';
funct_folder = '9resultTest';
path = fullfile(root_path, method_folder, funct_folder, 'time');
filename = 'ValdRuntime0.004Th5NNFold3.txt';
filepath = fullfile(path, filename);
data = importdata(filepath);
data = data{1};
% class(data)
strArr = strsplit(data);
inds = 4:5:length(strArr);
times = strArr(inds);
img_num = length(times);
time_cell = times(img_num);
time_char = time_cell{1};
time = time_char(1:length(time_char)-2);
time_ave = 1.0 * time*60/img_num
