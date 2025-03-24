function intergratesign(path)
    subpath1=strcat(path,'snpm+/');
    subpath2=strcat(path,'snpm-/');
    system(strcat('fslmaths "',subpath1,'snpm+_filtered.nii" -nan "',subpath1,'snpm+_filtered_nan.nii"'));
    system(strcat('fslmaths "',subpath2,'snpm-_filtered.nii" -nan "',subpath2,'snpm-_filtered_nan.nii"'));
    system(strcat('fslmaths "',subpath1,'snpm+_filtered_nan.nii" -sub "',subpath2,'snpm-_filtered_nan.nii" "',path,'snpm.nii"'));
end
