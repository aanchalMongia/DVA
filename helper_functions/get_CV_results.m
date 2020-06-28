function [ auc,aupr,XcROC,YcROC,XcPR,YcPR, T ] = get_CV_results(Y,n,cv_setting,predictionMethod  )
%GET_CV_RESULTS Summary of this function goes here
%   Detailed explanation goes here

global Sd Sv
            AUCs  = zeros(1,n);  AUPRs = zeros(1,n);   
            % loop over the n folds
            
            for i=1:n

               [test_ind,len]=get_test_ind( cv_setting,i,n, Y,Sd,Sv);

               y2 = Y;
               y2(test_ind) = 0;
               fprintf('*');
               st=tic;
               y3=alg_template(y2,predictionMethod,test_ind ,[]);
               endt= double(toc-st);
               [AUCs(i),XcROC,YcROC]  = calculate_auc (y3(test_ind),Y(test_ind));    
               [AUPRs(i),XcPR,YcPR, T] = calculate_aupr(y3(test_ind),Y(test_ind));
               %if (length(test_ind)==(round(len/n)-1)) %for 1st fold test_ind has 1 less element than the total elements in test_ind from other 9 folds..so copy 1 ele to X coordinate and Y coordinate
                   %disp('Ìn ifffff')
               %end
               
            end
            auc= mean(AUCs);    aupr= mean(AUPRs);
            
end

