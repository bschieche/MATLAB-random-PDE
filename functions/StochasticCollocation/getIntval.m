function intval = getIntval(nodes,intInfo,d,rule)

% newintval = getIntval(nodes,intInfo,d)
%
% function to calculate multidimensional integrals of Lagrange polynomials
%
% nodes: nodes for which the integrals are needed
% intInfo: information about 1D integrals
% d: dimension
%
% intval: vector of results as long as nodes
%
% Bettina Schieche, 2nd of August 2010
% intInfo is now a cell array of intInfos for several 1D quadrature nodes, rule is needed to specify which dimension belongs to which
    % rule, 2nd of February 2011

n = size(nodes,1); % number of nodes
intval = zeros(n,1); % initialize result
for l = 1:n % loop over nodes
    int = 1;
    for i = 1:d % loop over all directions
        currentIntInfo = intInfo{rule{end}(i)};
        [~,j] = tolismember(nodes(l,i),currentIntInfo(:,1),10^-8); % find index of current node in intInfo
		int = int*currentIntInfo(j,2); % multiply with 1D-integral value of dimension i
    end
    intval(l) = int; % save value in intval
end
