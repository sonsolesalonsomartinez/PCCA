%% Analysis Figure 5
%Tests whether Functional connectivity trajectories are not driven by
%simple activation patterns using cross-validated regularized ridge
%regression on each PC time series separately
addpath(genpath('../MATLAB/NetsPredict-master'))

nICs = 25; cond_all = {'rest','WM','motor'}; nCond = length(cond_all);nPCs = 5;

% Model parameters
family = 'gaussian';
    
figure;
for cond = cond_all

    if strcmp (cond,'rest')
        nFrames = 1200;  nSess = 4; nSub = 100;
    elseif strcmp (cond, 'WM')
        nFrames = 405;   nSess = 2; nSub = 99;
    elseif strcmp (cond,'motor')
         nFrames = 284;   nSess = 2; nSub = 100;
    end
    
    fileData = [cond '/tsData.mat'];  % timeseries data
    fileHmm0 = [cond '/hmm/hmmResults'];
    filePca0 = [cond '/pcca/pccaResults' ];
    filePca1 = [filePca0 '1_200hmmRuns']; % first solution of all pcca runs
    
    % Response variable    
    load(filePca1,'scores');
    % Predictors
    ttrial = nFrames*nSess;
    load(fileData,'X');
    % Interaction term
    str = '[';
    for i = 1:nICs-1
        for j = i+1:nICs
            str = [str, 'X(:,' num2str(i) ').*X(:,' num2str(j) '),']; %#ok<AGROW>
        end
    end
    str=[str, ']'];
    XX = eval(str);
    % Model parameters
    correlation_structure = repelem(1:nSub,1,ttrial)';

    subplot(1,nCond,find(ismember(cond_all,cond)));
    %% 3.1.order-1 IC time series
    nbars=4;
    a=0;
    Xin = [ones(length(X),1), X]; % Add column of 1's to include constant term in regression
    parameters.Method = 'glmnet';
    for pc = 1:nPCs
        Yin = scores(:,pc);
        stats= nets_predict5(Yin,Xin,family,parameters,correlation_structure,[],[]);
        bar(a,stats.cod);hold on
        a=a+nbars;
    end
    %% 3.2. order-2 interactions
    a=1;
    Xin = [ones(length(X),1), XX];
    parameters.Method = 'glmnet';
    for pc = 1:nPCs
        Yin = scores(:,pc);
        stats= nets_predict5(Yin,Xin,family,parameters,correlation_structure,[],[]);
        bar(a,stats.cod);hold on
        a=a+nbars;
    end
    
    %% 3.3. order-2 interactions deconfounded by amplitude
    a=2;
    Xin = [ones(length(X),1), XX];
    confounds = X;
    parameters.Method = 'Ridge';
    for pc = 1:nPCs
        Yin = scores(:,pc);
        stats= nets_predict5(Yin,Xin,family,parameters,correlation_structure,[],confounds);
        bar(a,stats.cod);hold on
        a=a+nbars;
    end
    ylabel('variance explained (R^2)');
    xlabel('top 5 principal components')
    title(cond)
end
