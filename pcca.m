 function [coeff, scores, e,index] = pcca(fileHmm0,n,K,nPCs, nICs)
        
    total_nRuns=1000;

    M = zeros(nICs,n*K);

    % take n sets of K states from the pool of 1000 hmm runs
    rng shuffle
    index = randperm(total_nRuns,n);

    % concatenate M (n*k) HMM state time series (Gamma)
    for rep = 1:n          
        r = index(rep);
        gamma = cell2mat(struct2cell(load([fileHmm0 num2str(r) '.mat'], 'Gamma')));
        c=size(find(M(1,:)),2);
        M(:,c+1:c+K)=gamma;
    end

    % Run PCA
    [coeff,scores,e] = pca(M);

    % Keep the top 20 PCs for memory issues
    coeff  = coeff(:,1:nPCs);
    scores = scores(:,1:nPCs); 
end