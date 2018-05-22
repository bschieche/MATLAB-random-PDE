function intInfo = getIntInfo(deltagrid,a,b,varargin)

% intInfo = getIntInfo(deltagrid,a,b)
%
% function to calculate all possible integrals of 1D-Lagrange polynomials
%
% deltagrid: information about nodes on all levels
% [a,b]: domain
% varagin: if varargin is given, another 1D node structure is given (n+1 new nodes instead of n-1)
%
% intInfo: 1. column = node, 2. column = corresponding value of integral
%
% Bettina Schieche, 2nd of August 2010

maxnodeno = length(deltagrid{end,1}); % number of possible nodes
intInfo = zeros(maxnodeno,2); % initialize intInfo
l = 1; % counter
for i = 1:size(deltagrid,1) % loop over all nodes
    % get index set J for nodes which are new on level i, i=1,i=2: special cases
	if nargin>3 || i<=2
		J = 1:2:2.^i-1;
    else
        J = 2:2:(2^(i-1)+1);
    end
    for j = J % loop over all new nodes on level i
        node = deltagrid{i,1}(j); % extract the node
        intInfo(l,1) = node; % save it in intInfo
        % compute the integral -> constant weight of 1/(b-a) = density uniform distribution
        %intInfo(l,2) = 1/(b-a)*quadl(@(x)lagrange(x,deltagrid{i,1},node),a,b);
        intInfo(l,2) = 1/(b-a)*deltagrid{i,2}(j);
        l = l+1; % increase counter
    end
end
