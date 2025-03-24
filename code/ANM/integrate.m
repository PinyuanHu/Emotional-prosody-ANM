function integrate(path,targetpath)
    sub=dir(path);
    sub(1:2)=[];
    for i=1:(size(sub,1)/4)
        fname1=strcat(path,sub(4*i-3).name);
        fname2=strcat(path,sub(4*i-2).name);
        fname3=strcat(path,sub(4*i-1).name);
        fname4=strcat(path,sub(4*i).name);
        Y1=spm_read_vols(spm_vol(fname1));
        Y2=spm_read_vols(spm_vol(fname2));
        Y3=spm_read_vols(spm_vol(fname3));
        Y4=spm_read_vols(spm_vol(fname4));
        Y=Y1+Y2+Y3+Y4;
        Y(isnan(Y))=0;
        system(strcat('cp "',fname1,'" "',targetpath,'"'));
        V=spm_vol(strcat(targetpath,sub(4*i-3).name));
        spm_write_vol(V,Y);
        disp(fname1);
    end
end