path='/data1/pyhu/HCP/results_1084_20240709_gender/';
sub=dir(fullfile(path,'/LR/affective/fc_whole_female/*T*'));%T/A
newpath='/data1/pyhu/HCP/results_1084_20240709_gender/affective/female/tmaps/';
maskpath_AICHA='/data1/pyhu/HCP/results_1084_20240709_gender/affective/female/AICHA_masked_tmaps/';
maskpath_WM='/data1/pyhu/HCP/results_1084_20240709_gender/affective/female/WM_masked_tmaps/';
for i = 1:size(sub,1)
    disp(i);
    system(strcat('fslmaths "',path,'/LR/affective/fc_whole_female/',sub(i).name,'" -add "',path,'/RL/affective/fc_whole_female/',sub(i).name,'" -div 2 "',newpath,sub(i).name,'"'));
    
    system(strcat('fslmaths "',newpath,sub(i).name,'.gz" -mul "/data/home/pyhu/data/mask/AICHA/AICHA_mask.nii.gz" "',maskpath_AICHA,sub(i).name,'.gz"'));
    system(strcat('fslmaths "',newpath,sub(i).name,'.gz" -mul "/data/home/pyhu/data/mask/152seg_sub_WM.nii.gz" "',maskpath_WM,sub(i).name,'.gz"'));
    system(strcat('gunzip "',maskpath_AICHA,sub(i).name,'.gz"'));
    system(strcat('gunzip "',maskpath_WM,sub(i).name,'.gz"'));
end