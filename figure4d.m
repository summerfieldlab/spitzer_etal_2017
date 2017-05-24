%% Illustration of RSA-based 'neurometric mapping', forward model
% Bernhard Spitzer, 2016

clear all;
kaps=[0.5  2   2];  % exponents kappa (k)
offs=[ 0   0  0.3]; % offset bias (b)

figure;
for ka=1:length(kaps)

    kap=kaps(ka);
    off=offs(ka); 

    %% 'neurometric'
    f=[-1:0.4:1];
    neurom=abs((off+f).^kap).*sign(off+f);
    subplot(3,length(kaps),(ka-1).*length(kaps)+1);
    plot(neurom,'Color',[0.5 0.5 0.5],'Linewidth',1.5); hold on; plot(zeros(1,6),'k--'); xlim([0.5 6.5])

    %% distance (12 by 12; 6 red, 6 green)
    p=[neurom neurom];
    plotDV=[p(1:6);p(7:12)]';
    dist=[];
    for i=1:12
        for j=1:12
            dist(i,j)=abs(p(i)-p(j));
        end
    end
    subplot(3,length(kaps),(ka-1).*length(kaps)+2);
    imagesc(dist); colormap('hot'); 
    set(gca,'xtick',[]); set(gca,'ytick',[]);

    %% response-mapped distance (dist x category red/green; 12 by 12)
    p=[-neurom neurom];
    distXcat=[];
    for i=1:12
        for j=1:12
            distXcat(i,j)=abs(p(i)-p(j));
        end
    end
    subplot(3,length(kaps),(ka-1).*length(kaps)+3);
    imagesc(distXcat); colormap('hot');
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    
end

        
        
        
        