function [G pred g] = psymodfun(b,Y,X,ML,nk,f,gnorm)
    
%% logfunX_ml_power
% psychometric model for sequential comparison tasks
% (c)Bernhard Spitzer, 2016

    % weighting function
    pownum=abs((b(2)+X(:,1:nk)).^b(3)).*sign(b(2)+X(:,1:nk)); % b(2): bias; b(3): kappa
    
    % g  normalization
    g=sum(abs(f+b(2)).^b(3)) / sum(abs(f)); 
    if gnorm
        pownum=pownum/g;
    end
   
    % sign-flip according to sample category;
    pflip=pownum.*X(:,nk+1:2*nk); 
    
    % leakage term
    leaker=(1-b(5)).^(nk-[1:nk]');
    pfin=pflip*leaker; % if b(5)=0, this corresponds to sum(pflip,2)

    DV=b(1)+pfin; % decision value
    pred=1./(1+exp(-DV./b(4))); % logit rule - b(4): sigma (noise)
    
    % ML fit (not used in simulation)
    if ML==1 
       G=sum(2*log(1./pred(Y==1))) + sum(2*log(1./(1-pred(Y==0))));
    else
       G=sum((Y-pred).^2); % ols
    end
end
