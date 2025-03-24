% 20240305 WHOLE BRAIN GENE
% clc
% clear
% close all;

% Gene Association Analysis
% define work dir
data_dir = '/data1/pyhu/GENE/general_compare/';

% Atlas
pathA = '/data/home/pyhu/data/mask/AICHA/AICHA.nii';
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
% 获取 Atlas 中的唯一标签（假设 Atlas 是整数标签，每个标签对应一个区域）
unique_labels = 1:384;
% unique_labels = unique(atlas_nii(:));
% unique_labels = unique_labels(2:end,:);
% % 判断数组中的奇数并保留,left
% unique_labels = unique_labels(mod(unique_labels, 2) ~= 0);
% % 判断数组中的偶数并保留,right
% unique_labels = unique_labels(mod(unique_labels, 2) == 0);
% 创建一个array，用于存储每个区域对应的 T Map 值
tmap = zeros(size(atlas_nii));
j = 1;
% 遍历每个区域，找到对应的 T Map 值
for i = 1:length(unique_labels)
    label = unique_labels(i);
    if label == roi(j,:)
        % 在 Atlas 文件中找到当前区域的索引
        mask = (atlas_nii == label);
        % 提取对应区域的 T Map 值
        tmap(mask) = PLS1_scores(j,:);   
        j = j + 1;
    end
end
% 保存nii,一定要修改fname,不然原始数据会被覆盖掉
VA.fname = [data_dir,'BN_PLS1_scores.nii'];
VA.dt = [16 0];
spm_write_vol(VA,tmap); %保存为nii格式

% correlation between PLS score and t-map
[corr_val,p_val] = corr(Y,PLS1_scores);

% Draw figures for correlations
% close all
% dotcolor = [2 48 74]/255;
% linecolor = [0 0 0];
% ylable1 = 'PLS1 scores';
% xlable1 = {'{\itz}-statistic of hemisphere-specific map'};
% [xData, yData] = prepareCurveData(Y, PLS1_scores);
% ft = fittype( 'poly1' );
% opts = fitoptions( ft );
% opts.Lower = [-Inf -Inf];
% opts.Upper = [Inf Inf];
% [fitresult, gof] = fit( xData, yData, ft, opts );
% h=plot( fitresult, xData, yData);
% set(h(1),'Marker','.','MarkerSize',6,'Color',dotcolor)
% set(h(2),'LineWidth',0.5,'Color',linecolor)
% hold on
% xFit = linspace(min(xData),max(xData),100);
% yPredict = predint(fitresult,xFit,0.95,'functional','off');
% fy = cat(2,yPredict(:,2)',flip(yPredict(:,1),1)')';
% fx = cat(2,xFit,flip(xFit',1)')';
% fill(fx,fy,[0.5 0.5 0.5],'EdgeAlpha',0,'FaceAlpha',0.3);
% hold off
% legend off
% ylabel(ylable1);
% xlabel(xlable1);
% set(gca,'LineWidth',0.5);
% set(gca,'FontName','Arial','FontSize',10);
% t1 = text(-2.5,0.25,{'{\itr} = 0.376';'{\itp} < 0.0001'},'FontName','Arial','FontSize',10);
% grid off
% box off
% print(gcf,[data_dir,'corr_PLS1_tmap.tif'],'-dtiff','-r1000')
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
% 显著性p值
p_single = zeros(1,dim);
for l = 1:dim
    p_single(l) = length(find(PCTVARrand(:,l)>=PCTVAR(2,l)))/n_per;
end

myStats = [PCTVAR; p_single];
% PCTVAR是自变量和因变量的解释方差
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
