function PLS_gene_syw_BN(pathA,data_dir,n)
    %default AICHA
    VA = spm_vol(pathA);
    atlas_nii = spm_read_vols(VA);
    
    %% load t-map,gene
    MRIdataPath = [data_dir,'merge.csv'];
    MRIMatric = readmatrix(MRIdataPath, 'OutputType', 'string'); 
    MRIdata = str2double(MRIMatric(:,2));
    roi = str2double(MRIMatric(:,1));
    
    genePath = [data_dir,'expression.csv'];
    opts = detectImportOptions(genePath);
    opts.VariableNamingRule = 'preserve';
    genecsv=readtable(genePath,opts);
    genes = genecsv.Properties.VariableNames;
    genes = genes(2:end);
    GENEdata = table2array(genecsv(:,2:end));
    
    %% PLS_calculation
    X = GENEdata;
    Y = zscore(MRIdata);
    dim = 10;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,dim,'CV',15);
    
    % align PLS components with desired direction
    R = corr(XS(:,n),MRIdata);
    if R < 0
        XS(:,n) = -1*XS(:,n);
    end
    
    % PLS1 ROIscores tmap
    PLS_scores = XS(:,n);
    unique_labels = 1:384;
    tmap = zeros(size(atlas_nii));
    j = 1;
    for i = 1:length(unique_labels)
        label = unique_labels(i);
        if label == roi(j,:)
            mask = (atlas_nii == label);
            tmap(mask) = PLS_scores(j,:);   
            j = j + 1;
        end
    end
    
    VA.fname = strcat(data_dir,'BN_PLS',num2str(n),'_scores.nii');
    VA.dt = [16 0];
    spm_write_vol(VA,tmap); 
    
    % correlation between PLS score and t-map
    [corr_val,p_val] = corr(Y,PLS_scores);
    csvwrite(strcat(data_dir,'PLS',num2str(n),'_ROIscores.csv'),XS(:,n));
    
    % initial PLS
    [PLSw,x] = sort(stats.W(:,n),'descend');
    
    if n==1
        % permutation test
        n_per = 10000;
        PCTVARrand = zeros(n_per,dim);
        Rsq = zeros(n_per,1);
        parpool(15);
        parfor j = 1:n_per
        %     disp(j);
            order = randperm(size(Y,1));
            Yp = Y(order,:);
            [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr] = plsregress(X,Yp,dim);
            PCTVARrand(j,:) = PCTVARr(2,:); 
        end
        p_single = zeros(1,dim);
        for l = 1:dim
            p_single(l) = length(find(PCTVARrand(:,l)>=PCTVAR(2,l)))/n_per;
        end
        
        myStats = [PCTVAR; p_single];
        csvwrite([data_dir,'PLS_stats.csv'],myStats);
    end
    
    %% calculate corrected weight
    
    n_boot = 10000;
    dim = 1;
    
    PLSweights = zeros(size(X,2),n_boot);
    geneindex = 1:size(X,2);
    
    parfor i = 1:n_boot
        myresample = randsample(size(X,1),size(X,1),1);
        res(i,:) = myresample; % store resampling out of interest
        Xr = X(myresample,:); % define X for resampled regions
        Yr = Y(myresample,:); % define Y for resampled regions
        [XLb,YLb,XSb,YSb,BETAb,PCTVARb,MSEb,statsb]=plsregress(Xr,Yr,dim); % perform PLS for resampled data
    
        temp = statsb.W(:,1);% extract PLS1 weights
        newW = temp(x); % order the newly obtained weights the same way as initial PLS 
        if corr(PLSw,newW) < 0 % the sign of PLS components is arbitrary - make sure this aligns between runs
            newW = -1*newW;
        end
    
        PLSweights(:,i) = newW;% store (ordered) weights from this bootstrap run  
        disp(i);
    end
    
    % standerd erros
    PLSsw = std(PLSweights');
    temp = PLSw./PLSsw';
    [Z,ind] = sort(temp,'descend');
    PLS = genes(ind);
    geneindex = geneindex(ind);
    fid = fopen(strcat(data_dir,'PLS',num2str(n),'_geneWeights.csv'),'w');
    for i = 1:length(genes)
      fprintf(fid,'%s, %d, %f, %f\n', PLS{i}, geneindex(i), Z(i), abs(Z(i)));
    end
    fclose(fid);
