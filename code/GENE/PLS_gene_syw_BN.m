function PLS_gene_syw_BN(pathA,data_dir)
% Gene Association Analysis
% Atlas
VA = spm_vol(pathA);
atlas_nii = spm_read_vols(VA);

%% load t-map,gene
MRIdataPath = [data_dir,'merge.csv'];
MRIMatric = readmatrix(MRIdataPath, 'OutputType', 'string'); 
MRIdata = str2double(MRIMatric(:,2));
roi = str2double(MRIMatric(:,1));

genePath = [data_dir,'AICHA_gene_all.csv'];
opts = detectImportOptions(genePath);
opts.VariableNamingRule = 'preserve';
genecsv=readtable(genePath,opts);
genes = genecsv.Properties.VariableNames;
genes = genes(2:end);
GENEdata = table2array(genecsv(:,2:end));

%% PLS_calculation
X = GENEdata;
Y = zscore(MRIdata);
dim = 15;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,dim,'CV',10);

% plot(1:15,cumsum(100*PCTVAR(2,:)),'-bo');
% xlabel('Number of PLS components');
% ylabel('Percent Variance Explained in y');

% align PLS components with desired direction
R1 = corr(XS(:,1),MRIdata);
if R1 < 0
    XS(:,1) = -1*XS(:,1);
end

% PLS1 ROIscores tmap
PLS1_scores = XS(:,1);
unique_labels = 1:384;

tmap = zeros(size(atlas_nii));
j = 1;
for i = 1:length(unique_labels)
    label = unique_labels(i);
    if label == roi(j,:)
        mask = (atlas_nii == label);
        tmap(mask) = PLS1_scores(j,:);   
        j = j + 1;
    end
end
VA.fname = [data_dir,'BN_PLS1_scores.nii'];
VA.dt = [16 0];
spm_write_vol(VA,tmap);

% correlation between PLS score and t-map
[corr_val,p_val] = corr(Y,PLS1_scores);
csvwrite([data_dir,'PLS1_ROIscores.csv'],XS(:,1));

% initial PLS
[PLS1w,x1] = sort(stats.W(:,1),'descend');

%% permutation test
n_per = 10000;
PCTVARrand = zeros(n_per,dim);
Rsq = zeros(n_per,1);
parpool(15);
parfor j = 1:n_per
    disp(j);
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

%% calculate corrected weight

n_boot = 10000;
dim = 1;

PLS1weights = zeros(size(X,2),n_boot);
geneindex = 1:size(X,2);

parfor i = 1:n_boot
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:) = myresample; % store resampling out of interest
    Xr = X(myresample,:); % define X for resampled regions
    Yr = Y(myresample,:); % define Y for resampled regions
    [XLb,YLb,XSb,YSb,BETAb,PCTVARb,MSEb,statsb]=plsregress(Xr,Yr,dim); % perform PLS for resampled data

    temp = statsb.W(:,1);% extract PLS1 weights
    newW = temp(x1); % order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW) < 0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW = -1*newW;
    end
    PLS1weights(:,i) = newW;% store (ordered) weights from this bootstrap run  
end

% standerd erros
PLS1sw = std(PLS1weights');
temp1 = PLS1w./PLS1sw';
[Z1,ind1] = sort(temp1,'descend');
PLS1 = genes(ind1);
geneindex1 = geneindex(ind1);
fid1 = fopen([data_dir,'PLS1_geneWeights.csv'],'w');
for i = 1:length(genes)
  fprintf(fid1,'%s, %d, %f, %f\n', PLS1{i}, geneindex1(i), Z1(i), abs(Z1(i)));
end
fclose(fid1);
