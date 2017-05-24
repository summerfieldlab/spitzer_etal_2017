%% RSA-based 'neurometric' mapping
%(c) Bernhard Spitzer, 2016
%% 
clear all; restoredefaultpath;
addpath('thirdparty');
datapath='rdmdat\';

selmod=1; % 1: visual; 2: auditory 

subjects={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24'}; %N=24

kaps=[0.4:0.02:4]; % exponents kappa (k)
offs=[-1:0.02:1];  % offset bias (b)
% this can take a while; lower resolution (increasing 0.02 above) will be faster, with similar mean results 

modlab={'vis';'aud'}; 
lookt=[0.2 0.6];
modmap=[];

for sub=1:length(subjects)
    
    disp(['loading data file subject ' subjects{sub}])
    actdat=load([datapath modlab{selmod} 'RDM' subjects{sub} '.mat']);
    
    D=actdat.Dset; % read settings
    epoT=D.epoT;
    subnum=str2num(subjects{sub});

    mahalRSA(:,:,:,sub)=actdat.RDM;
    looksmp=find(epoT>=lookt(1) & epoT<=lookt(2));
    subRDM=mean(actdat.RDM(:,:,looksmp),3);
      
    %% map over kappas and offsets
    for ka=1:length(kaps)
        kap=kaps(ka);
        for of=1:length(offs)
            off=offs(of);
            [mods orthmods] = rdm_neurom(kap,off,0);
            models=fieldnames(orthmods); 
            for m=1:2 % 1st: dist; 2nd:  dist x cat;
                actmod=squareform(orthmods.(models{m+3})); % +3: skip stimulus models
                orthmodmap(ka,of,m,sub)=rankCorr_Kendall_taua(actmod,squareform(subRDM));
            end
        end
    end
    
end
morthmodmap=mean(orthmodmap,4);
mmorthmodmp=mean(morthmodmap,3);

%% individual subjects's ~best fits (e.g., for corr with behavior)
suborthpara=[]; colsuborthpara=[];
for n=1:length(subjects)
    actmod=mean(orthmodmap(:,:,:,n),3); 
    [C,I] = max(actmod(:));
    [ka,of] = ind2sub(size(actmod),I);
    suborthpara(n,1)=kaps(ka);
    suborthpara(n,2)=offs(of);   
    colsuborthpara(n)=C;
end

%% dist/dist x cat (conjunction)
thr=0.001;
testdat=orthmodmap;
refpnt=orthmodmap(kaps==1,offs==0,:,:); % symmetric-linear reference
mrefpnt=mean(mean(refpnt));
for i=1:size(testdat,1)
    for j=1:size(testdat,2)
        for k=1:size(testdat,3)
            [p h(i,j,k)]=signrank(squeeze(testdat(i,j,k,:)),0,'alpha', thr);
        end
    end
end

conjmat=NaN(size(h(:,:,1)));
alpha=0.7
conjmask=ones(size(h(:,:,1))).*alpha;
for k=1:size(h,1)
    for o=1:size(h,2)
        if sum(h(k,o,:))==2 & mean(mean(testdat(k,o,:,:),3),4)>mrefpnt 
           conjmat(k,o)=mean(mean(testdat(k,o,:,:),3),4);
           conjmask(k,o)=1;
        end
    end
end

figure; 
h=imagesc(offs,kaps,mmorthmodmp); colorbar; xlabel('b'); ylabel('kappa (k)'); hold on;
set(h, 'AlphaData', conjmask);
plot(offs,ones(1,length(offs)),'--','Color',[0.5 0.5 0.5],'LineWidth',2);
plot([0 0],[kaps(1),kaps(end)],'--','Color',[0.5 0.5 0.5],'LineWidth',2); 
scatter(suborthpara(:,2),suborthpara(:,1),[],colsuborthpara,'filled','LineWidth',1,'MarkerEdgeColor','k');   
caxis([0 max(max(mmorthmodmp))+0.015]);

