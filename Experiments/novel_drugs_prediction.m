clear all
predictionMethod = 'grmf'; 

%----add dependencies to path----
addpath(genpath('../helper_functions'));


%read data
load('../data_processed/virus_drug_association.mat')
mat=mat'; %size of data matrix: #drugsx#vir
global Sd Sv
load('../data_processed/drug_sim_matrix.mat')
load('../data_processed/vir_sim_matrix.mat')
Y=mat; St=Sv;

%----parameters----
getParameters(predictionMethod,3)
    
       
       drugs_singleAssociation_indices=find(sum(Y,2)==1);
       other_drugs_indces=find(sum(Y,2)~=1);
       
       y2=Y;
       y2( drugs_singleAssociation_indices ,:) = 0;
       
       M = ones(size(Y)); 
       M( drugs_singleAssociation_indices,:)=0;
       test_ind=find(M==0);
       
       fprintf('*');

       y3=alg_template(y2,predictionMethod,test_ind ,[]);
        
       %-----------------------------------------------
       count=0;
       for i =1:length(drugs_singleAssociation_indices)
           drug_index=drugs_singleAssociation_indices(i);
           drug_pred=y3(drug_index,:);
           drug_actual=Y(drug_index,:);
           
           top1disease_index=find(drug_pred==max(drug_pred));
           if(drug_actual(top1disease_index)==1)
               count=count+1;
           end
                  
       end
       count
       no_singleAss_drugs=length(drugs_singleAssociation_indices)
       perc_correct_pred=count*100/length(drugs_singleAssociation_indices)
       [AUC,a,a]  = calculate_auc(y3(test_ind),Y(test_ind)); 
       [AUPR,a,a,a] = calculate_aupr(y3(test_ind),Y(test_ind));
       
       