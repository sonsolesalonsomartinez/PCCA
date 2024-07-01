% Analysis Figure 2
% 
load(filePca1, 'index');
nRuns = 200;

%% Pannel A
% Pearson correlation coefficient between states for each pair of HMM runs
% (aligned to coincide in their ordering), for K = 12 states and 200 HMM
% runs (i.e. 19900 pairs).
CC=zeros(nRuns*(nRuns-1)/2,K);
count=0;
for r1 = 1:nRuns-1
    gamma1=cell2mat(struct2cell(load([fileHmm0 num2str(index(r1)) '.mat'],'Gamma')));
    for r2=r1+1:nRuns
        count=count+1;
        gamma2=cell2mat(struct2cell(load([fileHmm0 num2str(index(r2)) '.mat'],'Gamma')));
        [S,assig, gamma2] = getGammaSimilarity(gamma1,gamma2);
        for k=1:K
            g1=gamma1(:,k);             
            g2=gamma2(:,k); 
            CC(count,k)=corr(g1,g2);
        end
    end
end
figure;violinplot(CC);ylabel('correlation coefficient');xlabel('sorted HMM states');

%% Pannel B
%Plot representing the percentage of individual (bars) and cumulative
%(line) explained variance by the top 20 PCs that resulted from applying
%PCA to the K = 12 states from 200 HMM runs (i.e., 12 × 200 = 2400 states).
load(filePca1,'e');
eb = (e/sum(e))*100;
ec=cumsum(e)/sum(e);
ec=ec*100;

figure;
yyaxis left
bar(eb,'edgecolor','none','FaceColor',0.8*[1 1 1]);hold on
ylabel('individual variance (%)');
yyaxis right
plot(ec(1:nPca),'.-','Color',[0.3647    0.1294    0.3804],'linewidth',1.5,'MarkerSize',9);box off;hold on
xlim([0.5 nPca+0.5]);set(gca,'XTick',5:5:nPca,'YColor',[0.3647    0.1294    0.3804]);
ylabel('cummulative (%)');
xlabel('top 20 principal components')

%% Pannel C
%Compute the Pearson correlation coefficient (y-axis) of the 20 PCs
%(x-axis) between each pair of 30 PCA runs, for several numbers of HMM runs
%(from 20 to 200)
nPCCApairs = (nPCCA*(nPCCA-1))/2;
cc_results = nan(length(nRuns2test)*nPCCApairs,nPCs);
cc_results_idx = nan(length(nRuns2test)*nPCCApairs,1);
count=0;
for n = nRuns2test
    scoresn=zeros(nPCCA,nICs,nPCs);  
    for j = 1:nPCCA
        load([filePca0 num2str(j) '_' num2str(n) 'hmmRuns.mat'],'scores');
        scoresn(j,:,:) = scores(:,1:nPCs);
    end
    % compute corelations
    for j1 = 1:nPCCA-1
        scores_j1 = squeeze(scoresn(j1,:,1:nPCs));
        for j2 = j1+1:nPCCA    
            count = count +1;
            scores_j2 = squeeze(scoresn(j2,:,1:nPCs));
            S = corr(scores_j1, scores_j2);
            cc_results_idx(count,1)=n;
            cc_results(count,:)=diag(S);
        end
    end  
end

figure;
mycolor = parula(length(nRuns2test));colormap(mycolor)
x = linspace(1,nPCs,nPCs);
pcount=0;
for n = nRuns2test
    pcount = pcount+1;
    ind = find(cc_results_idx(:,1)==n);
    mean_cc = mean(abs(cc_results(ind,1:nPCs)));
    err_cc = std(abs(cc_results(ind,1:nPCs)))/sqrt(length(ind));
    errorbar(x,mean_cc,err_cc,...
            '-o','MarkerSize',1,'MarkerFaceColor',mycolor(pcount,:),...
            'MarkerEdgeColor',mycolor(pcount,:),'Color',mycolor(pcount,:));hold on
end
c=colorbar();c.Ticks=0:0.09:1;c.TickLabels=nRuns2test;title(c,sprintf('nRuns'));
ylabel('Between-PCA-runs correlations');xlabel('principal components');
ylim([0 1])
%% Pannel D
%Compute the average state contribution (PCA coefficients) of each HMM run
%to the final PCCA model for the top 20 PCs.
load(filePca1,'coeff');
W = nan(nRuns,nPCs);
acc=0;
for r = 1:nRuns
    W(r,:) = sum(coeff(acc+(1:K),1:nPCs).^2);
    acc=acc+K;
end
figure;imagesc(W);ylabel('HMM run #');xlabel('principal components')
c=colorbar;title(c,'weights');