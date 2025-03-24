function ANMoverlapmap(path,thr)
    sub=dir(path);
    sub(1:2)=[];
    mkdir(strcat(path,'_nan'));
    mkdir(strcat(path,'_nan_bin+'));
    mkdir(strcat(path,'_nan_bin-'));
    mkdir(strcat(path,'_nan_bin'));
    mkdir(strcat(path,'_nan_bin_overlap'));
    for i=1:size(sub,1)
        system(strcat('fslmaths "',path,'/',sub(i).name,'" -nan "',path,'_nan/',sub(i).name,'"'));
        system(strcat('fslmaths "',path,'_nan/',sub(i).name,'" -thr +',num2str(thr),' -bin "',path,'_nan_bin+/',sub(i).name,'"'));
        system(strcat('fslmaths "',path,'_nan/',sub(i).name,'" -add 100 -uthr +',num2str(100-thr),' -bin "',path,'_nan_bin-/',sub(i).name,'"'));
        system(strcat('fslmaths "',path,'_nan_bin+/',sub(i).name,'" -sub "',path,'_nan_bin-/',sub(i).name,'" "',path,'_nan_bin/',sub(i).name,'"'));
    end
    system(strcat('cp "',path,'_nan_bin/',sub(1).name,'.gz" "',path,'_nan_bin_overlap/merge.nii.gz"'));
    for i=2:size(sub,1)
        system(strcat('fslmaths "',path,'_nan_bin_overlap/merge.nii.gz" -add "',path,'_nan_bin/',sub(i).name,'" "',path,'_nan_bin_overlap/merge.nii.gz"'));
        disp(i)
    end
end