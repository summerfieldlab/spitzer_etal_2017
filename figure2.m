clear all;
%% define experiment/simulation
ntrials=10000;
figrow=3; %1: Fig.2a; 2: Fig.2b upper; 3: Fig.2b lower;

if figrow==1 % main experiment
    nk=10; % k samples (per sequence) 
    ns=[1:6]; % range of number samples 
    ofs=[0.4375 0.2668]; % offset bias (visual, auditory) 
    plotes=[1.9474 1.9115]; % exponent kappa (visual, auditory)
    ylims=[-1.5 2.5; -1.5 2.5]; % for weighting function plot
    acculim=[0.55 0.85; 0.55 0.9]; % for weighting function plot
    estnoise=[1.7514  1.8529]; % noise (s) estimates (vis, aud) 
elseif figrow==2 % supporting exp, standard
    nk=8; 
    ns=[1:9]; 
    ofs=[0.2567 0.4055];  % (dots, digits)
    plotes=[1.552 1.641]; % (dots, digits)
    ylims=[-1.5 2; -1.5 2]; 
    acculim=[0.55 0.95; 0.55 0.85]; 
    estnoise=[1.3283  1.3283]; % common s estimate
elseif figrow==3 % supporting exp, standard
    nk=8;
    ns=[1:9]; 
    ofs=[0.01572 0.1038]; % (dots, digits) 
    plotes=[2.79 4.004];  % (dots, digits)
    ylims=[-1.5 2; -1.5 2]; 
    acculim=[0.55 1; 0.55 1]; 
    estnoise=[4.8013  4.8013]; % common s estimate
end

%% some specs
plotf=-1:0.05:1;
const=0;
lincol=[0 76 153]./256;
estcol=[0.5 0 0.5];
optcol=[237 125 49]./256;
nlincol=estcol;
nplotrows=2;
f=-1:2/ns(end-1):1;

%% simulate uniform random data
rng('shuffle');
for i=1:ntrials
    nums(i,:)=randsample(ns,nk,1);        % numbers
    catdat(i,:)=randsample([-1 1],nk,1);  % categories (e.g. red/green)
end
numdat=(nums-mean(ns))./max(ns-mean(ns));
Xdat=numdat.*catdat;
Y=sign(sum(Xdat,2))./2+0.5; % -> 0 1
catsum=sum(catdat,2);
nozeros=find(ismember(Y,[0 1]) & abs(catsum)<10); % excludes indifferent trials 
X=[numdat catdat];
X=X(nozeros,:);
Ytru=Y(nozeros); 

%% map parameter space
prenoises=[0.2:0.2:1.2 1.5 2.0 2.5 0.5:4]; % general
es=[0.4:0.1:6]; % exponents (kappa k)
nex=length(es);
plotnex=length(plotes);

figure;
% simulate accuracy
for o=1:length(ofs) 
    noises=[prenoises estnoise(o)]; % add estimated noise to plot
    noises=sort(noises);
    estnindex=find(noises==estnoise(o)); % and remember it's index
    nno=length(noises); 
    
    accu=[]; accucomp=[]; fnum=[];
    for e=1:length(es)
        for noise=1:length(noises)
            setb=[const ofs(o) es(e) noises(noise) 0]; 
            [G Ypred g]=psymodfun(setb,Ytru,X,0,nk,f,1);
            accu(noise,e)=1-mean(abs(Ytru-Ypred));
        end
    end
    
    % weighting function
    fnum(1,:)=(abs((ofs(o)+plotf).^plotes(o)).*sign(ofs(o)+plotf));
    fdisc(1,:)=(abs((ofs(o)+f).^plotes(o)).*sign(ofs(o)+f));
    
    subplot(nplotrows,4,1+(o-1).*2);
    colormap gray; mycols=flipud(colormap); 
    mycols=mycols(10:end,:);
    colstep=length(mycols)./plotnex;
    colvec=floor(1:colstep:length(mycols));
    mycols=mycols(colvec,:);
    
    % plot weighting function
    set(0,'DefaultLineLineWidth',1.5)
    plot(plotf,fnum(1,:),'Color',nlincol); hold on; 
    plot(f,fdisc(1,:),'o','MarkerSize',4,'MarkerFacecolor',nlincol,'Color',nlincol); hold on; 
    set(0,'DefaultLineLineWidth',1.5)
    ylim(ylims(o,:)); xlim([-1.1 1.1]);
    set(0,'DefaultLineLineWidth',3)
    plot(f,zeros(1,length(f)),'k--');     
    xlabel('X'); ylabel('decision value');
    
    
    % plot accuracy simulation
    subplot(nplotrows,4,2+(o-1).*2); title(['of=' num2str(ofs(o))]);
    set(0,'DefaultLineLineWidth',3)
    col=[0 0 0];
    for i=1:nno
        plotdat=accu;
        if i==estnindex
            plot(es,plotdat(i,:),'LineWidth',2,'Color',estcol); xlabel('kappa'); ylabel('accuracy'); hold on;
        else
            plot(es,plotdat(i,:), 'Color',col); xlabel('kappa'); ylabel('accuracy'); hold on;
        end
        col=col+(1./(nno+1));
    end
    set(0,'DefaultLineLineWidth',3)
    [C I]=max(plotdat');
    outsidemax=find(I==length(es));
    if length(outsidemax)>1
        C(outsidemax(2:end))=NaN;
    end
    plot(es(I),C,'--','Color',optcol); hold on; % plot max in orange
    set(gca,'XTick',[0:1:6])
    xlim([0 es(end)+0.5]);
    ylim(acculim(o,:));
    xlabel('kappa (k)'); ylabel('accuracy'); 
    
    %% add marker for (i) linear k=1 and (ii) fitted kappa
    plot([plotes(o) plotes(o)],acculim(o,:),'--','LineWidth',2,'Color',estcol);

end
