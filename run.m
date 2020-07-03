clear all
predictionMethod = 'gr1bmc_ppxa'%'gr1bmc_ppxa'%'mc 'dmf'%'mf' %'grmc_admm' %'grmf' 

% read virus-drug assocaitions
load('data_processed/virus_drug_association.mat')
mat=mat'; %size of data matrix: #drugsx#vir

global Sd Sv
load('data_processed/drug_sim_matrix.mat')
load('data_processed/vir_sim_matrix.mat')
Y=mat; 

%----add dependencies to path----
addpath(genpath('helper_functions'));

%----define parameters----
n = 10;% 'n' in "n-fold experiment"
global f_roc f_pr

%tic
for cv_setting=[ 1 2 3 ] 
  
getParameters(predictionMethod,cv_setting)
[auc,aupr,XcROC,YcROC,XcPR,YcPR, T ]=get_CV_results(Y,n,cv_setting,predictionMethod  );
   
auc
aupr

end
toc