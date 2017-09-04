function [clusters] = clustperm(dat, thrsh, clustthr, nits, nonparam)
% 1-dimensional cluster-based permutation test against ~ 0
% dat: conditions x samples x subjects; clusters: conditions x samples; 
% (c) Bernhard Spitzer, 2016
    
    clusters=NaN(size(dat,1),size(dat,2));
    for i=1:size(dat,1)
        
        disp(['cluster-based permutation test ' num2str(i) '/' num2str(size(dat,1))]);
        
        itclust=zeros(1,nits);
        for iter=1:nits
            itstats=zeros(1,size(dat,2));
            flipper=randsample([-1 1],size(dat,3),1)';
            for t=1:size(dat,2)
                if nonparam
                    [p,h,stats] = signrank(squeeze(dat(i,t,:)).*flipper);
                else
                    [h p stats] = ttest(squeeze(dat(i,t,:)).*flipper);
                end
                itstats(t)=p;
            end

            % find clusters (iterations)
            array=find(itstats<thrsh);
            temp = abs(diff(array));
            indx = find(temp > 1);
            indx = [0 indx length(array)];
            iclust=zeros(length(indx)-1,1);
            for ii = 1:(length(indx)-1)
               iclust(ii) = length(array((indx(ii)+1):indx(ii+1)));
            end 
            itclust(iter)=max(iclust);
        end
        itclust=sort(itclust,2,'descend');
        ctoff=clustthr.*nits;
        thrclust=itclust(ctoff);

        actstats=zeros(1,size(dat,2));
        for t=1:size(dat,2)
            if nonparam
                [p,h,stats] = signrank(squeeze(dat(i,t,:)));
            else
                [h p stats] = ttest(squeeze(dat(i,t,:)));
            end
            actstats(t)=p;
        end

        % find clusters (data)
        array=find(actstats<thrsh);
        temp = abs(diff(array));
        indx = find(temp > 1);
        indx = [0 indx length(array)];
        actclust={};
        for ii = 1:(length(indx)-1)
           actclust{ii} = array((indx(ii)+1):indx(ii+1));
           if length(actclust{ii})>thrclust
               clusters(i,actclust{ii})=i;
           end
        end 
        
    end     

end

