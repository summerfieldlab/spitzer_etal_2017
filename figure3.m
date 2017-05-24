%% ERP/CPP, leave-one-out time windows

clear all; 
addpath('thirdparty');

%%
selmod=2;    % 1: visual; 2: auditory
l1o=1;       % leave-one-out stats (1: yes, 0: no)
savepath='C:\Users\Berni\Dropbox\numcum\ms\NHB_rev1\eeg_check\';
thrsh=0.01;  % FDR-threshold in leave-one-out testing for significant time-windows

%%
subjects={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24'}; %N=24

%%
modlab={'vis';'aud'};
for sub=1:length(subjects)
    
    disp(['loading data file subject ' subjects{sub}])
    actdat=load(['erpdat\' modlab{selmod} 'ERP' subjects{sub} '.mat']);
    
    D=actdat.Dset; % read settings
    goodch(:,sub)=~D.badchannels; % bad channel indices
    subnum=str2num(subjects{sub});
    
    % average over reds and greens (for basic erp analysis only)
    for num=1:6
        erpdat(:,:,num,sub)=mean(actdat.ERPdat(:,:,[num num+6]),3);
    end

end

% average over ROI, excluding badchannels
selchans={'CP1', 'P1', 'POz', 'Pz', 'CPz', 'CP2', 'P2'};
epoT=D.epoT;
chanind=find(ismember(D.chanlabels,selchans));
subdat=[];
for sub=1:length(subjects)
    subchind=find(ismember(chanind,find(goodch(:,sub))));
    subchan=chanind(subchind);
    subdat(:,:,sub)=squeeze(mean(erpdat(subchan,:,:,sub),1));
end

% leave-one-out
subinds=1:length(subjects);
myclust=50; % in samples
l1outTwin=[]; bardat=[]; subTwins=[];
%figure; 
for actout=1:length(subjects)
    actin=find(subinds~=actout);
    if ~l1o
        actin=subinds;
    end 
    actindat=subdat(:,:,actin);
    Fcourse=[]; p=[];
    for t=1:length(epoT)
        testdat=permute(actindat(t,:,:),[3 2 1]);
        [p(t),tbl,stats] = kruskalwallis(testdat,[],'off');
        Fcourse(t)=tbl{2,5};
    end
    [p_fdr, p_masked] = fdr(p,thrsh);
    supfin=find(p_masked==1);
    p_clustmasked=zeros(size(p_masked));
    for k=2:length(supfin) %  cluster find
        if supfin(k)+myclust<length(epoT)
            if  p_masked(supfin(k))==1 & sum(p_masked(supfin(k):supfin(k)+myclust-1))==myclust 
                p_clustmasked(supfin(k):supfin(k)+myclust-1)=1;
            end
        end
    end
    %plot(epoT,p_clustmasked+actout*2); hold on; %  plot l1o windows
    tmpl1out=find(p_clustmasked==1);
    bardat(:,actout)=squeeze(mean(subdat(tmpl1out(1):tmpl1out(end),:,actout)));
    subTwins(:,actout)=[tmpl1out(1) tmpl1out(end)];
end
  
%% bar plot
figure;
colormap parula; mycols=colormap;
mycols=flipud(mycols);
mycols=mycols(10:end,:);
colstep=length(mycols)./6;
colvec=floor(1:colstep:length(mycols));
mycols=mycols(colvec,:);

mbardat=mean(bardat,2);
SEbar=std(bardat,[],2)./sqrt(sub);
subplot(2,2,3)
for i=1:length(mbardat)
    bar(i,mbardat(i),'Facecolor',mycols(i,:)); hold on
end
errorbar(1:6,mbardat,SEbar,'.k','Linewidth',1); ylabel('normalized amplitude');
set(gca,'xtick',[1:6]);


%% Whole-head ERP and topography plotting
% % Requires Fieldtrip toolbox 
%
% %average over subjects, excluding badchannels
% merpdat=[]; SEerpdat=[];
% for ch=1:64
%     goodchsubs=find(goodch(ch,:));
%     merpdat(ch,:,:)=mean(erpdat(ch,:,:,goodchsubs),4);
%     SEerpdat(ch,:,:)=std(erpdat(ch,:,:,goodchsubs),[],4)./sqrt(sub);
% end
% 
% dat=[]; cfg=[];
% pre.label=D.chanlabels(1:64)';
% pre.dimord='chan_time';
% pre.time=epoT;
% for k=1:size(merpdat,3)
%     lab=char(strcat('k', num2str(k)));
%     dat.(lab)=pre;
%     dat.(lab).avg=merpdat(:,:,k);
% end
% cfg.layout='ordered';
% lay = ft_prepare_layout(cfg, dat.k1);
% [X,Y]=getcoords(D.chanlabels(1:64)');
% lay.pos(1:64,:)=[X;Y]';
% scaler=0.6; %workaround for a nicely scaled topoplot
% lay.pos(1:length(D.chanlabels(1:64)),:)=[X;Y]'.*scaler;
% lay.width=lay.width*scaler;
% lay.height=lay.height*scaler;
% cfg.comment=[];
% cfg.colormap='jet'
% cfg.layout= lay;
% cfg.interactive='yes';
% cfg.style='straight';
% cfg.graphcolor=[1 1 0; 1 0.8 0.2; 1 0.6 0.4; 1 0.4 0.6; 1 0.2 0.8; 1 0 1];
% cfg.linewidth=2  
% figure;
% ft_multiplotER(cfg, dat.k1, dat.k2, dat.k3, dat.k4, dat.k5, dat.k6);

%% Plot channel ROI
plotdat=mean(subdat,3);
plotSE=std(subdat,[],3)./sqrt(sub);

subplot(2,2,1)
set(gca, 'ColorOrder', mycols, 'NextPlot', 'replacechildren');
ylimiter=[min(min(plotdat))*1.3 max(max(plotdat))*1.15];
ylimiter=[-0.5 1];
plot([0 0],ylimiter,'k--','Linewidth',1); hold on
x=epoT.*1000;
shadT=[epoT(round(mean(subTwins(1,:)))) epoT(round(mean(subTwins(2,:))))]; % average leave-one-out time window for shading
rectpos=[shadT(1).*1000 min(ylimiter) (shadT(2)-shadT(1))*1000 max(ylimiter)-min(ylimiter)];
rectangle('Position',rectpos,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[1 1 1],'Clipping','off'); 
plot(x,zeros(1,length(x)),'k--','Linewidth',1);
for i=1:size(plotdat,2)
    y=plotdat(:,i);
    dy=plotSE(:,i);
    fill([x';flipud(x')],[y-dy;flipud(y+dy)],mycols(i,:),'linestyle','none'); alpha(.5);
    line(x,y,'Color',mycols(i,:),'Linewidth',3); hold on
end
xlabel('time since sample onset (ms)'); xlim([-100 900]); set(gca,'layer','top')
ylabel('normalized amplitude'); ylim(ylimiter); 




