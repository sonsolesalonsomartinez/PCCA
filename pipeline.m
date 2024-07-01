
%
% This script requires the HMM-MAR toolbox (available at
% https://github.com/OHBA-analysis/HMM-MAR) and nets_predict5.m
% (availabel at https://github.com/vidaurre/NetsPredict)
%
% This script assumes that all the preprocessing has been done and that the
% IC time series are stored within the corresponding condition-folder:
% /rest,/WM and /motor in tsData.mat
%
% tsData.mat contains:
%   X, The IC time series of the fMRI data
%      concatenated across subjects (nICs by [nSub × nFrames x nSess])
%   T, the length of each session. For example for the resting-state
%      data T is a cell and T{subject} = [1200 1200 1200 1200] because
%      there are 1200 frames per run.
% 
% The pipeline follows the structure of the article:
% 1. Estimates FC trajectories using the Principal COmponents of
%    Connectivity Analysis(PCCA)
% 2. Tests the simmilarities of the spatial patterns of the estimated FC
%    trajectories
% 3. Tests the relationship between WM variables and FC trajectories at
%    different timescales.
% 4. Tests how much of the temporal covariation in FC can be predicted by
%    changes in the raw BOLD signal
%
% This pipeline must be adapted to your particular configuration of files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../MATLAB/HMM-MAR'));

% Setup filenames and data variables
TR=0.72; nICs=25; nCond=3; cond_all={'rest','WM','motor'}; 

cond = ''; %Run for each condition separately
if strcmp (cond,'rest')
    nFrames = 1200;  nSess = 4; nSub = 100;
elseif strcmp (cond, 'WM')
    nFrames = 405;   nSess = 2; nSub = 99;
elseif strcmp (cond,'motor')
     nFrames = 284;  nSess = 2; nSub = 100;
end
fileData = [cond '/tsData.mat'];      % timeseries data
fileHmm0 = [cond '/hmm/hmmResults'];  % HMM results
filePca0 = [cond '/pcca/pccaResults'];% PCCA results
filePca1 = [filePca0 '1_200hmmRuns']; % PCCA of run1 to be used for analysis

%% 1. Estimates FC trajectories using PCCA

% Load Time series data
load(fileData, 'X','T');

% HMM parameters
K = 12;
options = struct();
options.K = K; % number of states 
options.order = 0; % no autoregressive components
options.zeromean = 1; % do not model the mean
options.covtype = 'full'; % full covariance matrix
options.useParallel = 1;
options.cyc = 100;
options.standardise = 1;
options.verbose = 1;
options.inittype = 'HMM-MAR';
options.initcyc = 10;
options.initrep = 1;

% Run and store multiple HMMs
total_nRuns = 1000;
for r = 1:total_nRuns
    [hmm, Gamma] = hmmmar(X,T,options);
    save([fileHmm0 num2str(r) '.mat'], 'hmm', 'Gamma');
end

% Run and store several PCAs 
nRuns2test = [20:10:100 125 150 200];% no. HMM runs in each PCA run 
nPCCA = 30; % number of PCA runs
nPCs = 20;
for j = 1:nPCCA  
    for n = nRuns2test      
       [coeff, scores, e,index] = pcca(fileHmm0,n,K,nPCs,nICs);
       save([filePca0 num2str(j) '_' num2str(n) 'hmmRuns.mat'],'coeff','scores','e','index');
    end 
end

% Figure 1
open('analysis_figure_2.m');

%% 2.Tests the simmilarities of the spatial patterns of the estimated FC trajectories
open('analysis_figure_3.m');

%% 3.Tests the simmilarities of the spatial patterns of the estimated FC trajectories
open('analysis_figure_4.m');

%% 4. Tests how much of the temporal covariation in FC can be predicted by changes in the raw BOLD signal
open('analysis_figure_5.m');


