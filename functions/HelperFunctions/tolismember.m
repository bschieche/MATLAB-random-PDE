function [c,j] = tolismember(node,nodelist,tol)

% function [c,j] = tolismember(node,nodelist,tol)
%
% The function looks if nodelist contains node up to tolerance tol.
% 
% node = row vector determining a collocation point
% nodelist = matrix of nodes (columnwise)
% tol = tolerance when nodes should be treated as equal
%
% c = 1 if node is a member of nodelist up to a tolerance tol, c = 0 else
% j = index of node in nodelist which equals node
%
% Bettina Schieche, 2010

n = size(nodelist,1); % number of nodes

% find all nodes in nodelist with distance to node smaller than tol
I = find(sum(abs(bsxfun(@minus,nodelist,node)),2)<=tol);
if isempty(I)
    	c = 0;
		j = [];
else
    	c = 1;
    	j = I(1);
end
    
    
