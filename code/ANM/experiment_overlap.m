function experiment_overlap(seedpath,atlaspath,output)

    atlasv=spm_vol(atlaspath);
    atlasy=spm_read_vols(atlasv);
    sub=dir(seedpath);
    sub(1:2)=[];
    num_roi=max(max(max(atlasy)));
    a=zeros(num_roi,1);
    for i=1:size(sub,1)
        subpath=fullfile(sub(i).folder,sub(i).name);
        subv=spm_vol(subpath);
        suby=spm_read_vols(subv);
        for j=1:num_roi
            subroi=suby(atlasy==j);
            roiflag=(subroi~=0);
            if sum(roiflag)~=0
                a(j,1)=a(j,1)+1;
            end
        end
    end
    writematrix(a/size(sub,1),output);
    % b=zeros(size(atlasy));
    % for k=1:num_roi
    %     b(atlasy==k)=a(k,1)/size(sub,1);
    % end
    % atlasv.fname=output;
    % spm_write_vol(atlasv,b);
end
