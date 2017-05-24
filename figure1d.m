%% Figure 1d
clear all;

optimist=optimset('MaxFunEvals',50000,'MaxIter',50000,'Display','off');
subjects={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24'};

runs=[1:6];
inclmods=[1 2]; % modalities (1, visual; 2, auditory);
mods={'Visually';'Auditory'}; 
fignames={'visual';'auditory'}; 

nk=10;      % samples/stream
ns=1:6;     % number range
consterm=1; % include constant term (0: no, 1: yes) 
letleak=0;  % model leakage (0: no, 1: yes) if 1, model comparison will be against non-linear without leakage
plotonly=1; % plots precise (non-binary) model predictions 
figrows=1;  % adjust as needed
figcolumns=2;
kmax=10;   
smax=8;

butts=repmat([1 2],1,ceil(length(subjects)./2));
f=-1:2/ns(end-1):1;
fig=figure;
set(fig,'defaultAxesColorOrder',[0 0 0;0.4 0.4 0.4]);
if ~plotonly
    rng('shuffle'); % initialize random generator (binomial drawing below)
end
joinpsym=[]; joinbeta=[]; joinbetalin=[]; perf=[];

for modal=1:length(inclmods)
    actmod=inclmods(modal);
    nbeta=[]; nsubjmetr=[]; noptkappa=zeros(1,length(subjects));
    for n=1:length(subjects)
        vismat=[]; audmat=[]; 
        for run=runs %% load data
            blocknum=round(run/2);
            subnum=str2num(subjects{n});
            if ismember(subnum,[1:2:15 16:2:24]); % condition order (vis/aud) got flipped between 15/16
               modsel=mod(run,2)+1;
            elseif ismember(subnum,[2:2:14 17:2:23]); 
               modsel=~mod(run,2)+1;
            end 
            actfile=['behavdat\Numcum_' mods{modsel} '_Sub_' subjects{n} 'Block_' num2str(blocknum) '_' num2str(run) '.mat'];
            actfile;
            load(actfile);

            resp=all_digits(:,11);
            resp(~ismember(resp,[10 6]))=999; % so will be excluded below
            resp=round(resp/6)-1; % [6 10] -> [0 1]
            if mod(subnum,2)
                resp=-resp+1; % according to left/right button press assignment 
            end
            
            if modsel==1
                vismat=[vismat; all_digits(:,1:nk) resp]; %11 is resp
            elseif modsel==2
                audmat=[audmat; all_digits(:,1:nk) resp];
            end
        end % loop data load per subject
        
        if actmod==1 & numel(actmod)==1;
            data=vismat;
        elseif actmod==2 & numel(actmod)==1;
            data=audmat;
        elseif actmod==[1 2]
            data=[audmat;vismat];
        end
        numdat=(abs(data(:,1:nk))-3.5)./2.5;    
        catdat=sign(data(:,1:nk));
        catsum=sum(catdat,2);
        Xdat=numdat.*catdat;
        Y=data(:,11);
        X=[numdat catdat Xdat];

        Ytru=sign(sum(Xdat,2))./2+0.5;
        nozeros=find(ismember(Y,[0 1]) & abs(catsum)<10);
       
        X=X(nozeros,:);
        Y=Y(nozeros);
        Ytru=Ytru(nozeros);
        trueq=(Ytru==0.5); % exclude mean(red)=mean(green) from calculation of percentage correct responses 
        if sum(unique(Y)~=0 & unique(Y)~=1)>0
            error('check response vectors');
        end
        tricount(n)=length(Y);
        perf(modal,n)=mean(Ytru(~trueq)==Y(~trueq));

        %% estimate
        offStart=0; kappaStart=1; noiseStart=1; leakStart=0;
        b0=[0  offStart kappaStart noiseStart leakStart];
        lb=[-inf  -1  0.1   0  0]; 
        ub=[inf   1  kmax   smax  letleak];
        if ~consterm
            lb(1)=0; ub(1)=0;
        end
        beta=fmincon(@(b) psymodfun(b,Y,X,1,nk,f,0),b0,[],[],[],[],lb,ub,[],optimist); % non-linear model
        nbeta(n,:)=beta;
        [Gopt(modal,n) pred g]=psymodfun(beta,Y,X,1,nk,f,0);
        nbeta(n,4)=nbeta(n,4)/g; % rescale s/g
        
        if ~letleak
            lb(3)=1; ub(3)=1; % to compare with linear model
        elseif letleak
            lb(5)=0; ub(5)=0; % to compare with non-linear model without leakage 
        end
        betalin=fmincon(@(b) psymodfun(b,Y,X,1,nk,f,0),b0,[],[],[],[],lb,ub,[],optimist);
        [Goptlin(modal,n) predlin g]=psymodfun(betalin,Y,X,1,nk,f,0);
        nbetalin(n,:)=betalin;
        nbetalin(n,4)=nbetalin(n,4)/g; % rescale s/g

        %% psychometric functions
        subjresp=repmat(Y,1,nk);

        if ~plotonly
            predresp=binornd(1,pred);
            linresp=binornd(1,predlin);
        else
            predresp=pred; 
            linresp=predlin; 
        end
        predresp=repmat(predresp,1,nk);
        linresp=repmat(linresp,1,nk);

        numbers=X(:,1:nk).*2.5+3.5;
        predresp=reshape(predresp',[prod(size(numbers))],1);
        subjresp=reshape(subjresp',[prod(size(numbers))],1);
        linresp=reshape(linresp',[prod(size(numbers))],1);
        numbers=reshape(numbers',[prod(size(numbers))],1);
        cats=reshape(X(:,nk+1:2*nk)',[prod(size(numbers))],1);
          
        for num=1:6
            predmetr(num,n)=mean([predresp(numbers==num & cats==1); -(predresp(numbers==num & cats==-1))+1]); 
            subjmetr(num,n)=mean([subjresp(numbers==num & cats==1); -(subjresp(numbers==num & cats==-1))+1]); 
            linmetr(num,n)=mean([linresp(numbers==num & cats==1); -(linresp(numbers==num & cats==-1))+1]); 
        end
        
    end
    
    join4SE=cat(3,subjmetr,predmetr,linmetr);
    SE=squeeze(std(join4SE,[],2))./sqrt(n);
    plotmean=squeeze(mean(join4SE,2));

    subplot(figrows,figcolumns,modal); % Psychometric
    ylimiter=[0.4 0.7];
    yyaxis right
    ONEcol=[0 76 153]./256; % Linear
    errorbar(plotmean(:,3),SE(:,3),'--','Linewidth',1.5,'Color',ONEcol); hold on
    datcol=[0 0 0]; % data
    errorbar(plotmean(:,1),SE(:,1),'Linewidth',1.5,'Color',datcol);% hold on
    plot(plotmean(:,1),SE(:,1),'Linewidth',2,'Color',datcol); hold on
    MODcol=[112 48 160]./256; % non-linear (purple)
    plot(plotmean(:,2),'.','Linewidth',1,'MarkerSize',25,'Color',MODcol); %hold on;
    xlabel('number'); ylabel('choice probability') 
    title(fignames{actmod}); xlim([0.5 6.5]); 
    ylim(ylimiter);
    yyaxis left
    plot(0.*ones(1,6),'k--','Linewidth',1.5); ylabel('decision weight'); ylim(ylimiter-0.5);
    if modal==1
        set(gca,'xticklabel',{'1','2','3','4','5','6'}); 
        xlabel('digit');
    elseif modal==2
        set(gca,'xticklabel',{'one','two','three','four','five','six'},'XTickLabelRotation',0);
        xlabel('number word');
    end
    
    joinbeta=cat(3,joinbeta,nbeta(:,2:4)); % collect parameter estimates (visual, auditory)
    
end

meanbetas=squeeze(mean(joinbeta,1)) % columns visual, auditory; rows b, k, s;

