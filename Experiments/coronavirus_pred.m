clear all

predictionMethod = 'grmc_admm'

%----add dependencies to path----
addpath(genpath('../helper_functions'));

%----read data---- %
load('../data_processed/virus_drug_association_withcv.mat')
mat=mat'; %size of data matrix: #drugsx#vir
global Sd Sv
load('../data_processed/drug_sim_matrix_withcv.mat')
load('../data_processed/vir_sim_matrix_withcv.mat')
Y=mat; St=Sv;

getParameters(predictionMethod,2)
 %-----------------------------------------
k=5
           vir_ind=find(strcmp(vi_names,'SARS-CoV-2'));

           y2=Y;
           y2( :,vir_ind ) = 0;

           M = ones(size(Y)); 
           M( :,vir_ind)=0;
           test_ind=find(M==0);

           fprintf('*');

           y3=alg_template(y2,predictionMethod,test_ind ,[]);
          
           drugs_pred=y3(:,vir_ind);
           [sortedValues,sortIndex] = sort(drugs_pred,'descend');
           ind = sortIndex(1:k);
       
           top5drugNames=dr_names(ind)'
           sortedValues(1:5)'
