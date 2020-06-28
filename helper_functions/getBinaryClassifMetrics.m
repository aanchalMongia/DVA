function [precision,recall,Acc,F1] = getBinaryClassifMetrics(yHaT,yval)
%GETBINARYCLASSIFMETRICS Summary of this function goes here
%   Detailed explanation goes here

        tp = sum((yHaT == 1) & (yval == 1));
        fp = sum((yHaT == 1) & (yval == 0));
        fn = sum((yHaT == 0) & (yval == 1));
        tn = sum((yHaT == 0) & (yval == 0));
        
        precision = tp / (tp + fp);
        recall = tp / (tp + fn);
        F1 = (2 * precision * recall) / (precision + recall);
        Acc=(tp+tn)/(tp + tn + fp+ fn);

        
        
        %[c,cm,ind,per] = confusion(Y(test_ind),y3(test_ind)); 
      %{
       precision0 =  cm(2,2) / sum(cm(2,:));
        recall0 =  per(2,3) ;
       
        precision = mean(@(cm) diag(cm)./sum(cm,2));
        recall = mean(@(cm) diag(cm)./sum(cm,1)');
        f1Scores = @(cm) 2*(precision(cm).*recall(cm))./(precision(cm)+recall(cm));
        meanF1 = @(cm) mean(f1Scores(cm));
       %}
end

