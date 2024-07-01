%% Analysis Figure 3
%Tests the simmilarities of the spatial patterns of the estimated FC
%trajectories
Isupdiag=find(triu(ones(nICs,nICs),1));

% Get and transform the covariance matrices of the 20 PCs of each
% behavioural condition (rest, WM and motor) into correlation matrices and
% then applied the Fisher z-transformation on the off-diagonal elements of
% the matrices;
for cond = cond_all

    if strcmp (cond,'rest')
        nFrames = 1200;  nSess = 4; nSub = 100;
    elseif strcmp (cond, 'WM')
        nFrames = 405;   nSess = 2; nSub = 99;
    elseif strcmp (cond,'motor')
         nFrames = 284;   nSess = 2; nSub = 100;
    end

    fileHmm0 = [cond '/hmm/hmmResults'];
    filePca0 = [cond '/pcca/pccaResults' ];
    filePca1 = [filePca0 '1_200hmmRuns']; % first solution of all pcca runs

    c_cond=nan(nPCs,length(Isupdiag)); % cond1: REST
    C = getCovariancePC(filePca1,nPCs,nICs,fileHmm0);
    for i=1:nPCs
        Ci = squeeze(C(i,:,:));
        % transform cov matrices into corr matrices
        Ci = corrcov(Ci,0);
        %Fisher transformation
        c_cond(i,:) = atanh(Ci(Isupdiag)); 
    end
    eval(['c_' cond '= c_cond'])

end

% Compared the vectorised elements of the off diagonal to one another
% using Pearson correlation
C_wi_Rest = zeros(nPCs,nPCs);
C_wi_WM = zeros(nPCs,nPCs);
C_wi_Motor = zeros(nPCs,nPCs);
C_btw_RestWM = zeros(nPCs,nPCs);
C_btw_RestMotor = zeros(nPCs,nPCs);
C_btw_WMMotor = zeros(nPCs,nPCs);
for pc1 = 1:nPCs  
    for pc2 = 1:nPCs
        C_wi_Rest(pc1,pc2) = min(min(corrcoef(c_rest(pc1,:),c_rest(pc2,:))));
        C_wi_WM(pc1,pc2) = min(min(corrcoef(c_WM(pc1,:),c_WM(pc2,:))));
        C_wi_Motor(pc1,pc2) = min(min(corrcoef(c_motor(pc1,:),c_motor(pc2,:))));
       
        C_btw_RestWM(pc1,pc2) = min(min(corrcoef(c_rest(pc1,:),c_WM(pc2,:))));
        C_btw_RestMotor(pc1,pc2) = min(min(corrcoef(c_rest(pc1,:),c_motor(pc2,:))));
        C_btw_WMMotor(pc1,pc2) = min(min(corrcoef(c_WM(pc1,:),c_motor(pc2,:))));
    end
end

%% Pannel A
M = [C_wi_Rest  C_btw_RestWM  C_btw_RestMotor;
     C_btw_RestWM' C_wi_WM  C_btw_WMMotor;
     C_btw_RestMotor' C_btw_WMMotor' C_wi_Motor];

% Store correlation values in the lower triangle of the matrix
Isubd=find(triu(ones(length(M))));
CClow=M;
CClow(Isubd)=0;
% Store partial correlation values in the upper triangle of the matrix
CCup=-inv(M);
CCup=(CCup ./ repmat(sqrt(abs(diag(CCup))),1,nPCs*nCond)) ./ repmat(sqrt(abs(diag(CCup)))',nPCs*nCond,1);
Isupd=find(tril(ones(length(CCup))));
CCup(Isupd)=0;
% Plot matrix
figure;
CC2 = CClow+CCup;
imagesc(CC2,[0.2 1]);
set(gca,'XTick',0.5:5:nPCs*nCond,'XTickLabel',...
    '', 'YTick',0.5:5:nPCs*nCond,'YTickLabel','');
axis square;grid on;colorbar();hold on;

%% Pannel B
% Test the statistical significance of the averaged similarity of the
% spatial patterns (of rows with respect to columns) for all pairs of
% within (rest, WM, motor) and between (rest-WM, rest-motor, WM-motor)
% conditions. Based on permutation testing (10000 permutations).
cmp_labels={'rest','WM','motor','rest-WM','rest-motor','WM-motor'};
nCmp = length(cmp_labels);
Isub=find(tril(ones(nPCs),-1));
nPerm=10000;

pval=2*ones(nCmp,nCmp);
side=zeros(nCmp,nCmp);

i=1;j=2;[p,~,s] = permtestdiff(C_wi_Rest(Isub),C_wi_WM(Isub),nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 
i=1;j=3;[p,~,s] = permtestdiff(C_wi_Rest(Isub),C_wi_Motor(Isub),nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 
i=1;j=4;[p,~,s] = permtestdiff(C_wi_Rest(Isub),C_btw_RestWM,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end
i=1;j=5;[p,~,s] = permtestdiff(C_wi_Rest(Isub),C_btw_RestMotor,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 
i=1;j=6;[p,~,s] = permtestdiff(C_wi_Rest(Isub),C_btw_WMMotor,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end

i=2;j=3;[p,~,s] = permtestdiff(C_wi_WM(Isub),C_wi_Motor(Isub),nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 
i=2;j=4;[p,~,s] = permtestdiff(C_wi_WM(Isub),C_btw_RestWM,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 
i=2;j=5;[p,~,s] = permtestdiff(C_wi_WM(Isub),C_btw_RestMotor,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 
i=2;j=6;[p,~,s] = permtestdiff(C_wi_WM(Isub),C_btw_WMMotor,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 

i=3;j=4;[p,~,s] = permtestdiff(C_wi_Motor(Isub),C_btw_RestWM,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 
i=3;j=5;[p,~,s] = permtestdiff(C_wi_Motor(Isub),C_btw_RestMotor,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 
i=3;j=6;[p,~,s] = permtestdiff(C_wi_Motor(Isub),C_btw_WMMotor,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 

i=4;j=5;[p,~,s] = permtestdiff(C_btw_RestWM,C_btw_RestMotor,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end
i=4;j=6;[p,~,s] = permtestdiff(C_btw_RestWM,C_btw_WMMotor,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end

i=5;j=6;[p,~,s] = permtestdiff(C_btw_RestMotor,C_btw_WMMotor,nPerm);pval(i,j)=p;pval(j,i)=p;if p<0.05,side(i,j)=s;side(j,i)=-1*s;end 

% FDR corrections
Isub=find(tril(ones(nCmp),-1));
pval_vector=pval(Isub);
[~,~,~, pvals_bh_fig3]=fdr_bh(pval_vector,0.05);  


