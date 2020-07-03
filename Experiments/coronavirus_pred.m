clear all

predictionMethod = 'grmf'

%----add dependencies to path----
addpath(genpath('../helper_functions'));

%----read data---- %
load('../data_processed/withCorIsolates/virus_drug_association_withcv.mat')
mat=mat'; %size of data matrix: #drugsx#vir
global Sd Sv
load('../data_processed/withCorIsolates/drug_sim_matrix_withcv.mat')
load('../data_processed/withCorIsolates/vir_sim_matrix_withcv.mat')
Y=mat; St=Sv;

col0=find(strcmp(vi_names,'SARS-CoV-2')); Y(:,col0)=0;
col1=find(strcmp(vi_names,'SARS-CoV-2: feb'));  Y(:,col1)=0;
col2=find(strcmp(vi_names,'SARS-CoV-2: april'));  Y(:,col2)=0;
col3=find(strcmp(vi_names,'SARS-CoV-2: june'));  Y(:,col3)=0;

getParameters(predictionMethod,2)
 %-----------------------------------------
k=10
           vir_ind=find(strcmp(vi_names,'SARS-CoV-2'));
           %vir_ind=find(strcmp(vi_names,'SARS-CoV-2: feb'));
           %vir_ind=find(strcmp(vi_names,'SARS-CoV-2: april'));
           %vir_ind=find(strcmp(vi_names,'SARS-CoV-2: june'));

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
           sortedValues(1:5)';
