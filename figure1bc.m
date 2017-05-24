%% sim_powerf_accuracy
% Simulate psychometric observer and plot predicted model accuracy under varying 
% noise (sigma), offset bias (b), and kappa (k) parameters
% (c)Bernhard Spitzer, 2016

clear all; 
addpath(genpath('C:\Users\Berni\Dropbox\numcum\ms\NHB_rev1\upload\'));

%% set experiment/simulation

nk=10; % k samples (per sequence)
ns=[1:6]; % range of number samples 
eqcat=0; % if 1: equally many reds and greens


%% parameter space
noises=[0.2:0.2:1.2 1.5:0.5:4]; % range of noise levels
ofs=[0 0.2]; % offset biases (b)
es=[0.4:0.1:4]; % range of kappas (exponents k)


%% simulate uniform random data
ntrials=1000;
cats=[-ones(1,nk/2) ones(1,nk/2)]; % if eqcat
rng('shuffle');
for i=1:ntrials
    nums(i,:)=randsample(ns,nk,1); % numbers
    if eqcat
        catdat(i,:)=cats(randperm(nk)); % categories (red/green)
    else
        catdat(i,:)=randsample([-1 1],nk,1);
    end
end

numdat=(nums-mean(ns))./max(ns-mean(ns)); % rescale to -1...1
Xdat=numdat.*catdat; % sign-flip according to category (red/green)
Y=sign(sum(Xdat,2))./2+0.5; % ground truth (0: red, 1: green)

catsum=sum(catdat,2);
nozeros=find(ismember(Y,[0 1]) & abs(catsum)<10); % excludes indifferent trials 
X=[numdat catdat];
X=X(nozeros,:);
Ytru=Y(nozeros); 

%% some stuff for plotting
f=-1:2/ns(end-1):1;
plotf=-1:0.05:1;
plotes=[0.3:0.1:1 1.2:0.2:2 2.4:0.4:4]; % kappa levels (for line plots)
plotnes=length(plotes);
ONEind=find(plotes==1);
nno=length(noises);
nex=length(es);
const=0;
ylims=[-1.4 1.4; -1.5 3];
figrows=2;  %
figcolumns=3;

figure;

for o=1:length(ofs) 
    
    %% compute and plot weighting curve illustration
    for e=1:plotnes
        fnum(e,:)=(abs((ofs(o)+plotf).^plotes(e)).*sign(ofs(o)+plotf));
        fdisc(e,:)=(abs((ofs(o)+f).^plotes(e)).*sign(ofs(o)+f));
    end
    set(0,'DefaultLineLineWidth',2)
    subplot(figrows,figcolumns,o);
    colormap gray; mycols=flipud(colormap); 
    mycols=mycols(10:end,:);
    colstep=length(mycols)./plotnes;
    colvec=floor(1:colstep:length(mycols));
    mycols=mycols(colvec,:); 
    for i=1:plotnes
        if i~=ONEind;
            set(0,'DefaultLineLineWidth',2)
            plot(plotf,fnum(i,:),'Color',mycols(i,:)); hold on; 
            plot(f,fdisc(i,:),'o','MarkerSize',5,'MarkerFacecolor',mycols(i,:),'Color',mycols(i,:)); hold on; 
        end
    end
    set(0,'DefaultLineLineWidth',4)
    ONEcol=[0 76 153]./256;
    plot(plotf,fnum(ONEind,:),'Color',ONEcol); hold on; 
    plot(f,fdisc(ONEind,:),'o','Color',ONEcol,'MarkerSize',5,'MarkerFacecolor',ONEcol); hold on; 
    ylim(ylims(o,:)); xlim([-1.1 1.1]);
    set(0,'DefaultLineLineWidth',2)
    plot(f,zeros(1,length(f)),'k--');     
    xlabel('X'); ylabel('decision value');
    set(gca,'XTick',[-1 -0.5 0 0.5 1]);
    
    
    %% simulate and plot accuracy
    accu=[]; accucomp=[]; fnum=[];
    for e=1:length(es)
        for noise=1:length(noises)
            setb=[const ofs(o) es(e) noises(noise) 0]; 
            [G Ypred g]=psymodfun(setb,Ytru,X,0,nk,f,1);
            accu(noise,e)=1-mean(abs(Ytru-Ypred));
        end
    end
    
    if o==1
        subplot(figrows,figcolumns,o+2); title(['of=' num2str(ofs(o))]);
        set(0,'DefaultLineLineWidth',3)
        col=[0 0 0];
        for i=1:nno
            plotdat=accu;
            plot(es,plotdat(i,:), 'Color',col);  hold on;
            col=col+(1./(nno+1));
        end
        set(0,'DefaultLineLineWidth',3)
        [C I]=max(plotdat');
        outsidemax=find(I==length(es));
        if length(outsidemax)>1
            C(outsidemax(2:end))=NaN;
        end
        plot(es(I),C,'--','Color',[237 125 49]./256); hold on; % orange
        set(0,'DefaultLineLineWidth',1)
        xlim([0 es(end)+0.5]);
        if o==1
           ylim([0.55 1]);
        else
            ylim([0.6 0.95]);
        end
        xlabel('kappa (k)'); ylabel('accuracy'); 
        set(gca,'XTick',[0:1:max(es)])
    end

end
