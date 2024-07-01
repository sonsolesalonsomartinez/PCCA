function C = getCovariancePC(filePca,nPCs,nICs,basenameResults,K)
    
    load(filePca,'scores','index');
    nrep = length(index);
    nstates=K*nrep;
    scores=scores(:,1:nPCs);
    
    % compute covariance matrix for 
    C_pca = zeros(nstates,nPCs,nICs,nICs);
    counter=0;
    for i = 1:nrep
        ind = index(i);
        hmm = cell2mat(struct2cell(load(sprintf('%s%d.mat',basenameResults,ind), 'hmm')));
        gamma = cell2mat(struct2cell(load(sprintf('%s%d.mat',basenameResults,ind), 'Gamma')));
        for k = 1:K
            counter = counter+1;
            [C1, ~] = getFuncConn(hmm,k,1);
            for pc = 1:nPCs 
                beta = (gamma(:,k)' * (scores(:,pc).^2))/ sum(gamma(:,k));
                C_pca(counter,pc,:,:) = C1*beta;
            end   
        end 
    end
    C = squeeze(sum(C_pca));


end
