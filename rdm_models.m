function [mods orthmods] = rdm_models
%(c) Bernhard Spitzer, 2016
 
    %%  digit
    pc=ones(6,6);
    for i=1:6
        for j=1:6
            if i==j  % remove diagonal (6x6)
                pc(i,j)=0;
            end
        end
    end
    c=[pc pc; pc pc];
    mods.dig=c;

    %% category
    c=[zeros(6,6) ones(6,6); ones(6,6) zeros(6,6)];
    for i=1:12
        for j=1:12
            if i==j  % remove diagonal
                c(i,j)=0;
            end
        end
    end
    mods.cat=c;

    %%  parity
    v1=repmat([0 1],1,6);
    v2=repmat([1 0],1,6);
    c=repmat([v1;v2],6,1);
    mods.eve=c;

    %% numdist
    pc=zeros(6,6);
    for i=1:6
        for j=1:6
            pc(i,j)=abs(i-j)./5;
        end
    end
    c=[pc pc; pc pc];
    mods.numd=c;

    %%  numdist x category
    c=[pc fliplr(pc); fliplr(pc) pc];
    mods.nXc=c;

    %% orthogonalize models
    models=fieldnames(mods); % attention: check order above
    vecs=[];
    for m=1:length(models)
        vecs(:,m)=squareform(mods.(models{m}));
        vecs(:,m)=vecs(:,m)-mean(vecs(:,m)); % mean center
    end

    for m=1:length(models) % recursive orthogonalization
        order4orth=[find((1:length(models))~=m) m];
        tmpvecs=spm_orth(vecs(:,order4orth));
        orthmods.(models{m})=squareform(tmpvecs(:,end));
    end

%%    plot model RDMs
%     figure; colormap('hot');
%     for m=1:length(models)
%         subplot(2,length(models),m); 
%         imagesc(mods.(models{m}),[0 1]); title(models{m});
%         subplot(2,length(models),m+length(models));
%         imagesc(orthmods.(models{m}),[-1 1]); title([models{m} '_{orth}']);
%     end
    
end

