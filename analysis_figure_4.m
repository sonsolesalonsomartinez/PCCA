% Analysis Figure 4
% CCA to correlate block wise indices of PCCA (mean
% and variance over blocks of the WM paradigm) with task related measures
% (WM-load, RT and accuracy). This is done after empirical mode
% decomposition allowing us to see which frequencies of the PCAs
% are best related to the behavioral signal.

cond='WM'; nFrames = 405;   nSess = 2; nSub = 99; nBlocks = 8; TR = 0.72;
filePca0 = [cond '/pcca/pccaResults' ];
filePca1 = [filePca0 '1_200hmmRuns'];% first solution of all pcca runs

nPCs=5;
fs=1/TR;
TRms=TR*1000;
T1=1:nFrames;
T2=nFrames+1:nFrames*nSess;

load(filePca1,'scores');

% Calculate IMFs and their frequencies
imf={};
InstFreq=[];
for pc=1:nPCs
    
    [imf_pc, ~] = emd(scores(:,pc));
    nmodes=size(imf_pc,2);
    
    %% Pannel A
    % Plot IMFs from an illustrative session obtained from PC1’s time
    % series
    if pc==1
        strips(imf_pc(1:nFrames,:))
        set(gca,'YTickLabels',1:nmodes,'XTick',[0, nFrames-1],...
            'XTickLabels', round(TR*([1,nFrames])));
        ylabel('IMFs');xlabel('time (seconds)');
        title('IMFs of PC1 for a sample session')
    end
    
    for m = 1:size(imf_pc,2)  
        % Caluclate the Hilbert spectrum of the IFS
        [~, ~, ~,InstFreq(pc,m,:)] = hht(imf_pc(:,m),fs); %#ok<SAGROW>
        % Store all IMFs for each PC time series
        imf{m}{pc}=imf_pc(:,m); %#ok<SAGROW>
    end
end

%% Calculate CCA
% Load vars
vars  = dlmread('vars_blocks.txt',' ');% WM vars for each block and subject
twins = dlmread('twins.txt',' ');      % family structure

% Yin: WM RT and ACC 
Acc =  vars(:,1); % accuracy per subject and session and block
RT =   vars(:,2); % RT per subject and session and block
Load = vars(:,3); % indicates whether each block is 0-Back or 2-Back
T = [vars(:,4),vars(:,5)];% indicates the first and last time point period for each block
ps = repelem(1:nSub,nSess*nBlocks); % idx for each subject
Yin=[Load,Acc, RT]; % size(nSub*nSess*nBlocks,3)

% Xin: IMF mean and energy
nmodes=3;
Xin=nan(nSub*nSess*nBlocks,nPCs*2,nmodes);
acc=0;
for i = 1:length(Yin)
    tt = T(i,1):T(i,2);
    for m = 1:nmodes
        imf_mode = cell2mat(imf{m});
        Xin(i,1:nPCs,m) = mean(imf_mode(tt,:));
        Xin(i,nPCs+1:nPCs*2,m) = var(imf_mode(tt,:));
    end
end

%% WITHIN SUBJECT EFFECT
% canonical correlation permutations within Subjects
Y=Yin;
permPm=[];Rm=[];Am=[];Bm=[]; % mean
permPv=[];Rv=[];Av=[];Bv=[]; % energy
for m=1:nmodes
    
    % mean
    X=squeeze(Xin(:,1:nPCs,m));
    [permPm(m),~,Rm(m,:),~,~,~,Am(m,:,:),Bm(m,:,:)] = ...
        permtestcca0(X,Y,[],[],[],[],[],ps);%within
    
    % energy
    X=squeeze(Xin(:,nPCs+1:nPCs*2,m));
    [permPv(m),~,Rv(m,:),~,~,~,Av(m,:,:),Bv(m,:,:)] = ...
        permtestcca0(X,Y,[],[],[],[],[],ps);%within
end
% FDR
[~,~,~, pvals_bh]=fdr_bh([permPm,permPv],0.05);

% Figure 4C
figure;
bar(1./([pvals_bh(1:nmodes)',pvals_bh(nmodes+1:end)']),'EdgeColor','none','BarWidth',1);hold on;
yline(1/0.05,'--','p-FDR<0.05');hold on
yline(1/0.01,'--','p-FDR<0.01');hold on
set(gca,'YScale','log','box','off',...
    'YTick',[1,10,100,1000,10000],'YTickLabel',1./[1,10,100,1000,10000],...
    'YGrid','on','YMinorGrid','off','XTickLabel',compose('IMF%d',1:nmodes)); 
ylabel('p-FDR');
title('within-subjects effect')
legend({'mean','energy'},'box','off');

%% BETWEEN SUBJECT EFFECT
% canonical correlation permutations between subjects 
block2 = find(Yin(:,1) == 1);%(only high cognitive load blocks) 
Y=Yin(block2,2:3);

permPm=[];Rm=[];Am=[];Bm=[];% mean
permPv=[];Rv=[];Av=[];Bv=[];% energy
for m=1:nmodes
    % mean
    X=squeeze(Xin(block2,1:nPCs,m));
    [permPm(m),~,Rm(m,:),~,~,~,Am(m,:,:),Bm(m,:,:)] = ...
        permtestcca0(X,Y,[],[],[],twins,[],ps(block2));
    
    %energy
    X=squeeze(Xin(block2,nPCs+1:nPCs*2,m));
    [permPv(m),~,Rv(m,:),~,~,~,Av(m,:,:),Bv(m,:,:)] = ...
        permtestcca0(X,Y,[],[],[],twins,[],ps(block2));    
end
% FDR
[~,~,~, pvals_bh]=fdr_bh([permPm,permPv],0.05);

% Figure 4D
figure;
bar(1./([pvals_bh(1:nmodes)',pvals_bh(nmodes+1:end)']),'EdgeColor','none','BarWidth',1);hold on;
yline(1/0.05,'--','p-FDR<0.05');hold on
yline(1/0.01,'--','p-FDR<0.01');hold on
set(gca,'YScale','log','box','off',...
    'YTick',[1,10,100,1000,10000],'YTickLabel',1./[1,10,100,1000,10000],...
    'YGrid','on','YMinorGrid','off','XTickLabel',compose('IMF%d',1:nmodes)); 
ylabel('p-FDR');
title('between-subjects effect')
legend({'mean','energy'},'box','off');









