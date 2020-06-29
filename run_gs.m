clear all

%----add dependencies to path----
addpath(genpath('helper_functions'));

% read virus-drug assocaitions
load('data_processed/virus_drug_association.mat')
mat=mat'; %size of data matrix: #drugsx#vir
global Sd Sv
load('data_processed/drug_sim_matrix.mat')
load('data_processed/vir_sim_matrix.mat')
 
n=10;
%read data
Y=mat; 

global k alpha  k1 k2  pp lamda mu1 mu2    
global num_iter p lambda_l lambda_d lambda_t    


for cv_setting=[ 1 ] % 1 ds for 1 viral class

    %MC : No parameter tuning needed
    %{
    % MF
    predictionMethod = 'mf'
    bestauc = -Inf;
    l=[];
    for k=[35 30 25 20 15 10 5 ]
        for alpha=[0.01 0.1 0.5 1 1.5 2 ]
            
            [auc, aupr]=get_CV_results(Y,n,cv_setting,predictionMethod  )
            
            l=[l [k alpha auc aupr]'];
            save(['GS/AUCnAUPR_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'l')
            if bestauc < auc
                bestauc = auc;
                bestcomb = [k alpha];
               save(['GS/BestResult_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'bestcomb','aupr','bestauc')
            end
        end
    end
    %clear k alpha %clearvars instead Y Sd St
    
    
    
    %DMF
    predictionMethod = 'dmf'
    bestauc = -Inf;
    l=[];
    for alpha=[0.001 0.01 0.1 0.5 1 2   ]
      for k1=[35 30 25 20 ]
        for k2=[20 15 10 5 ]
           
            [auc, aupr]=get_CV_results(Y,n,cv_setting,predictionMethod  )
            
            l=[l [alpha k1 k2 auc aupr]'];
            save(['GS/AUCnAUPR_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'l')
            if bestauc < auc
                bestauc = auc;
                bestcomb = [alpha k1 k2];
               save(['GS/BestResult_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'bestcomb','aupr','bestauc')
            end
        end
      end
    end
    %clear k1 k2 %clearvars instead Y Sd St
    
 %}
    
    
    %GRMF
    
    predictionMethod = 'grmf'
    %global num_iter p k lambda_l lambda_d lambda_t
    bestauc = -Inf;
    l=[];
    num_iter = 2;
    k = 50;
     for p=2:7
         for lambda_l=2.^(-5:-1)
            for lambda_d=[0.005 0.01  0.05  0.1  0.15  0.2  0.25  0.3]
                for lambda_t=[0.005 0.01  0.05  0.1  0.15  0.2  0.25  0.3]
                      
                        [auc, aupr]=get_CV_results(Y,n,cv_setting,predictionMethod  )

                        l=[l [p lambda_l lambda_d lambda_t  auc aupr]'];
                        save(['GS/AUCnAUPR_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'l')
                        if bestauc < auc
                            bestauc = auc;
                            bestcomb = [p lambda_l lambda_d lambda_t];
                           save(['GS/BestResult_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'bestcomb','aupr','bestauc')
                        end    
                        
                    end
                end
            end
     end
     %}
     %{
    %GRDMF
    predictionMethod = 'dgrdmf'
    %global pp alpha k1 k2 mu1 mu2
    bestauc = -Inf;
    l=[];
    for pp=[2 5 8 ]
       for alpha=[ 0.01 0.1 1 2 ]
          for k1=[35 30 25 20 ]
            for k2=[20 15 10 5 ]
                for mu1=[0.001 0.01 0.1 1 10 100 ]
                    for mu2=[0.001 0.01 0.1 1 10 100]

            % loop over the n folds
            
             [auc, aupr]=get_CV_results(Y,n,cv_setting,predictionMethod  )
            
            l=[l [pp alpha k1 k2 mu1 mu2 auc aupr]'];
            save(['GS/AUCnAUPR_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'l')
            if bestauc < auc
                bestauc = auc;
                bestcomb = [pp alpha k1 k2 mu1 mu2];
               save(['GS/BestResult_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'bestcomb','aupr','bestauc')
            end
            
            end
           end
          end
         end
       end
    end
    %clear pp alpha k1 k2 mu1 mu2
    
    
%}
end

     
     
 for cv_setting=[1]        
    %GRMC_PPXA
    predictionMethod = 'gr1bmc_ppxa'
    %global pp lamda mu1 mu2 
    bestauc =0.8085;% -Inf;
    l=[];
    for pp=[2]
       for lamda=[0.05 2]%[0.1 0.5 0.05  ]
         for mu1= [0.01 0.05  0.1 1  ]
           for mu2=[0.1 0.5 1 ]
               
            [auc, aupr]=get_CV_results(Y,n,cv_setting,predictionMethod  )
            
            l=[l [pp lamda mu1 mu2  auc aupr]'];
            save(['GS/AUCnAUPR_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'l')
            if bestauc < auc
                bestauc = auc;
                bestcomb = [pp lamda mu1 mu2 ];
               save(['GS/BestResult_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'bestcomb','aupr','bestauc')
            end
            
           end
         end
       end
    end
    %clear pp lamda mu1 mu2  
    
 end 
 
 for cv_setting=[1  ] 
    %GRMC_ADMM
    predictionMethod = 'grmc_admm'
    %global  pp lamda mu1 mu2 
    bestauc = -Inf;
    l=[];
    for pp=[2 ]
       for lamda=[ 0.01 0.05 0.1 ]
         for mu1=[0.05 0.1 1 ]
           for mu2=[0.005 0.01 0.1 ]
            
            [auc, aupr]=get_CV_results(Y,n,cv_setting,predictionMethod  )
            
            l=[l [pp lamda mu1 mu2  auc aupr]'];
            save(['GS/AUCnAUPR_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'l')
            if bestauc < auc
                bestauc = auc;
                bestcomb = [pp lamda mu1 mu2 ];
               save(['GS/BestResult_method=' predictionMethod '_cv=' num2str(cv_setting) '.mat'],'bestcomb','aupr','bestauc')
            end
            
           end
         end
       end
    end
    %clear pp lamda mu1 mu2 

end
 %}
