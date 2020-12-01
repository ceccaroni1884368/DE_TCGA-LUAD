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
% ?!?!?
dataNDe = log2(dataNDe);
dataCDe = log2(dataCDe);

%% Parameter 
Z_th = 5; %  Z threshold

%% Correlation analysis and  co-expression network figure
%
[N, pv_N] = corr(dataNDe'); %compute Pearson's linear correlation coefficient
N = N-triu(tril(N));
[C, pv_C] = corr(dataCDe'); %compute Pearson's linear correlation coefficient
C = C-triu(tril(C));

fisherZN = atanh(N); %compute fisher z-transformation 
fisherZC = atanh(C); %compute fisher z-transformation 

[nn, nN] = size(dataNDe);
[nc, nC] = size(dataCDe);
ZScore = (fisherZN - fisherZC)/sqrt(1/(nN-3)+1/(nC-3));  %compute zscore

ZScore(abs(ZScore)<Z_th) = 0; %selecting only absolute values above the threshold

adjCN = logical(ZScore);
G_CN = graph(adjCN);
figure;h1 = plot(G_CN,'Layout','force','UseGravity',true);
title({'Differentially expressed genes';'Differentially Co-expression network'})


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
gene_symDe(HubsFC>0)

%% Subgraph
%  2    {'ZMYND10'}
% 15    {'NGEF'   }
% 67    {'NEK2'   }
% 84    {'CDHR3'  }
% 96    {'PEBP4'  }
%120    {'PLPP2'  }
%171    {'SFTPC'  }
%175    {'SFTPB'  }
%183    {'BUB1'   }
%209    {'MAP7D2' }
CentrMeasures = 'degree';
node = 209;
sG = subgraph(G, [node; neighbors(G,node)]);
index = centrality(sG,CentrMeasures);
Y = prctile(index,95);
Hubs_ind = find(index>Y);
HubsFC = mFCDe(Hubs_ind);
figure; h3 = plot(sG,'Layout','force','UseGravity',true);
highlight(h3,Hubs_ind,'NodeColor','r','Marker','h','MarkerSize',4)
labelnode(h3,[1:numnodes(sG)], gene_symDe([node; neighbors(G,node)]))

if find(HubsFC<0)
    highlight(h3,Hubs_ind((HubsFC<0)),'NodeColor','cyan','Marker','h','MarkerSize',4)
end

title(gene_symDe(node))

gene_symDe(HubsFC<0)
gene_symDe(HubsFC>0)