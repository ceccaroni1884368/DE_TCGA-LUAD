clc
clearvars
close all

%% Data loading
load('../TCGA-LUAD/data/data_DE.mat')
% dataNDe: Normal tissue data of Differentially expressed genes (DEGs)
% dataCDe: Cancer tissue data of DEGs
% gene_IDDe: gene IDs of DEGs
% gene_symDe: gene symbols of DEGs
% fdrt: FDR max (parameter used to define DEGs)
% fct: minimun fold change (parameter used to define DEGs)
% mFCDe: average Fold Chance of DEGs among the patients

%% Apply log2 to dataNDe and dataCDe
% ?!?
dataNDe = log2(dataNDe);
dataCDe = log2(dataCDe);

%% Parameter 
PCC_th = .7; % Pearson Correlation Coefficient (PCC) threshold

%% Correlation analysis and  co-expression network figure
% Normal data analysis
[CN, pv_CN] = corr(dataNDe'); %compute Pearson's linear correlation coefficient
CN(abs(CN)<PCC_th) = 0; %selecting only absolute values above the threshold

CN = CN-triu(tril(CN));
adjCN = logical(CN);
G_CN = graph(adjCN);
figure;h1 = plot(G_CN,'Layout','force','UseGravity',true);
title({'Differentially expressed genes';'Co-expression network in normal cells'})

%% Degree distribution
degree_CN = centrality(G_CN,'degree');
degree_CN(degree_CN==0) = [];
figure;histogram(degree_CN)    
title('Degree distribution')
xlabel('k');ylabel('Occurrences')
figure;histogram(degree_CN,'Normalization','probability')
title('Degree distribution')
xlabel('k');ylabel('P(k)')

%% Centrality Measures
%verify the presence of isolated nodes and elimination of them
K = sum(adjCN,2);
adjCN((K==0),:) = []; adjCN(:,(K==0)) = [];
G = graph(adjCN);
gene_IDDe((K==0)) = [];
gene_symDe((K==0)) = [];
mFCDe((K==0)) = [];
clear adjCN K

CentrMeasures = {'degree','closeness','betweenness','eigenvector'};

for cm=1:length(CentrMeasures)
    clear index Y
    index = centrality(G,CentrMeasures{cm});
    Y = prctile(index,95);
    figure; h2 = plot(G,'Layout','force','UseGravity',true);
    highlight(h2,(index>Y),'NodeColor','r','Marker','h','MarkerSize',4)
    title(CentrMeasures{cm})
    Hubs.(CentrMeasures{cm}) = gene_symDe(index>Y);
end

%% Compare hubs sets related to the use of the different centrality measures 
C = combnk(CentrMeasures,2);

for r=1:length(C)
    if r==1
        gene_intersect = intersect(Hubs.(C{r,1}),Hubs.(C{r,2}));
        if isempty(gene_intersect)
            break
        end
    else
        gene_intersect = intersect(gene_intersect,intersect(Hubs.(C{r,1}),Hubs.(C{r,2})));
        if isempty(gene_intersect)
            break
        end
    end
end

%% Only considering degree index, are the hubs up or down regulated?
CentrMeasures = 'degree';
index = centrality(G,CentrMeasures);
Y = prctile(index,95);
Hubs_ind = find(index>Y);
HubsFC = mFCDe(Hubs_ind);
figure; h3 = plot(G,'Layout','force','UseGravity',true);
highlight(h3,Hubs_ind,'NodeColor','r','Marker','h','MarkerSize',4)
if find(HubsFC<0)
    highlight(h3,Hubs_ind((HubsFC<0)),'NodeColor','cyan','Marker','h','MarkerSize',4)
end

title(CentrMeasures)

gene_symDe(HubsFC<0)