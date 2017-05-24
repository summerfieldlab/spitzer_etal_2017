%% RDM correlation time-courses / Grand mean RDM / MDS  
%  Bernhard Spitzer, 2016
%% 
clear all; restoredefaultpath;
addpath('thirdparty');
datapath='rdmdat\';

selmod=1; % 1: visual; 2: auditory 

subjects={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24'}; %N=24

thrsh=0.01;      % cluster-forming threshold
clustthr=0.01;   % cluster threshold
nits=1000;       % permutation test iterations
lookt=[0.2 0.6]; % display time window (grand-mean RDM and MDS) 

tic; mrsamat=[]; modcont=[];
modlab={'vis';'aud'};
for sub=1:length(subjects)
    
    disp(['processing EEG-RDM file subject ' subjects{sub}])
    actdat=load([datapath modlab{selmod} 'RDM' subjects{sub} '.mat']);
    
    D=actdat.Dset; % read settings
    epoT=D.epoT;
    subnum=str2num(subjects{sub});

    mahalRSA(:,:,:,sub)=actdat.RDM;
      
    %% obtain model RDMs
    [mods orthmods] = rdm_models;
    models=fieldnames(orthmods); % attention: check  order
    for m=1:length(models)
        actmod=squareform(orthmods.(models{m}));
        for t=1:length(epoT)
            modcont(m,t,sub)=rankCorr_Kendall_taua(actmod,squareform(actdat.RDM(:,:,t)));
        end
    end
    
end

%% permutaton test
rng('shuffle'); 
clusters=clustperm(modcont,thrsh,clustthr,nits,1);

%% plot RDM correlation time courses
plotdat=mean(modcont,3);
plotSE=std(modcont,[],3)./sqrt(sub);
plotdat=permute(plotdat,[2 1]); % just to fit plotting tool below
plotSE=permute(plotSE,[2 1]);
figure; 
colormap lines; mycols=colormap;
set(gca, 'ColorOrder', mycols, 'NextPlot', 'replacechildren');
ylimiter=[-0.1 0.28];
x=epoT.*1000;
for i=1:size(plotdat,2)
    y=plotdat(:,i);
    line(x,y,'Color',mycols(i,:),'Linewidth',3); hold on
end
legend(fieldnames(orthmods));
for i=1:size(plotdat,2)
    plot(x,0.006.*clusters(i,:)+ylimiter(1),'Linewidth',5,'Color',mycols(i,:)); % stats marker
    y=plotdat(:,i);
    dy=plotSE(:,i);
    fill([x';flipud(x')],[y-dy;flipud(y+dy)],mycols(i,:),'linestyle','none'); alpha(.5); hold on
end
plot([0 0],ylimiter,'k--','Linewidth',1); 
plot(x,zeros(1,length(x)),'k--','Linewidth',1);
xlabel('time since sample onset (ms)'); periT=[epoT(1) epoT(end)]; xlim(periT.*1000); 
ylabel('mean Tau'); ylim(ylimiter);

%% plot grand mean EEG RDM
looksmp=find(epoT>=lookt(1) & epoT<=lookt(2));
mRSA=mean(mahalRSA,4);
figure; 
clims=[0.64 0.84];
for wind=1:size(lookt,1)
    looksmp=find(epoT>=lookt(wind,1) & epoT<=lookt(wind,2));
    subplot(1,size(lookt,1),wind);imagesc(mean(mRSA(:,:,looksmp),3),clims); colormap('hot'); colorbar; 
end

%% MDS Figure
figure;
[mds e] = cmdscale(mean(mRSA(:,:,looksmp),3));
txoff=0.003;
plot3(mds(1:6,1),mds(1:6,2),mds(1:6,3),'.-r'); hold on;
text(mds(1:6,1)+txoff,mds(1:6,2)+txoff,mds(1:6,3)+txoff,cellstr(num2str([1:6]')),'Color','r','FontSize',16);
plot3(mds(7:12,1),mds(7:12,2),mds(7:12,3),'.-g')
text(mds(7:12,1)+txoff,mds(7:12,2)+txoff,mds(7:12,3)+txoff,cellstr(num2str([1:6]')),'Color','g','FontSize',16);
%%%%
