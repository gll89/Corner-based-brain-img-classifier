clear all, close all, clc;

prmts.window = 3;    %control the size of Gaussion window 
% rootPath = 'E:\BMCImgs\CTImgs';
prmts.root_path = '/home/labuser/LinlinGao/BMC/Image_revised';
prmts.categry_folder = 'CTImgs';
prmts.exts = 'bmp';   %'tiff';
prmts.ou_folder = '9result';
prmts.rslt_txt = 'cornerResultRecord'
prmts.subfolders = ['normal'; 'abnrml']

for th = 0.01 : 0.01 : 0.09     
    %For CT images, th belongs to [0.01, 0.09];  For MRI images, th belongs to [0.03, 0.12]
    prmts.th = th
    functCornerDetectMTIFinal(prmts);
end