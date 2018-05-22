function err = localError(newsurplus,newintval,surplus,intval,qoi,varargin)

% function to calculate local error indicator
%
% newsurplus = hierarchical surpluses for new nodes
%
% Bettina Schieche, 27th July 2010
% added more local errors, 2nd of August 2010

I = qoi;
addval = newintval'*newsurplus{1}(:,I);
meanval = intval'*surplus{1}(:,I);
tempErr = abs(addval)./abs(meanval); % maximal total contribution to several quantities of interest of integral type
err = max(tempErr);


