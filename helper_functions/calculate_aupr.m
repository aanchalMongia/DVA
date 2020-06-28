function [AUPR, X, Y,T]=calculate_aupr(predictionScores,labels)
%calculate_aupr calculates the area under the precision-recall curve
%
% INPUT:
%  predictionScores:    prediction scores
%  labels:              actual labels
%
% OUTPUT
%  AUPR:                area under the precision-recall curve
%
% Borrowed from code of:
% Twan van Laarhoven, Sander B. Nabuurs, Elena Marchiori,
% (2011) Gaussian interaction profile kernels for predicting drug�target interaction
% http://cs.ru.nl/~tvanlaarhoven/drugtarget2011/
 
global f_pr
%    figure(f_pr);
%    [X,Y,T,AUC] = perfcurve(labels,predictionScores,1 ,'XCrit','reca','YCrit','prec');
X=[]; Y=[]; T=[];
%    plot(X,Y,'Color',[0.85 0.85 0.85]);
%    hold on;
    
	if nargin > 1
		[~,i] = sort(predictionScores,'descend');
		labels = labels(i);
	end
	
	%for i=1:n
	%	if targets(i)
	%		goods = goods + 1
	%       aupr = aupr + good/(good+bad);
	%	else
	%		bad = bad + 1;
	%	end
	%end
	cumsums = cumsum(labels)./reshape(1:numel(labels),size(labels));
	AUPR = sum(cumsums(~~labels));
	pos = sum(labels);
	AUPR = AUPR / pos;
end