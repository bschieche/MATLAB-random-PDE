function basis = node2basis(nodelist,node,levellist,level,tol,flag)

% basis = node2basis(nodelist,node,levellist,level,tol)
%
% The function searches in all dimensions for nodes which contribute to the basis function of a certain node
%
% nodelist = list of nodes
% node = row vector determining one node
% levellist = information about the level of nodelist, the same size as nodelist
% level = information about the level of node
% tol = tolerance when nodes should be treated as equal
%
% basis = cell array of vectors with node indices for all dimensions
%
% Bettina Schieche, 21th of July 2010
%
% gets now a flag that reads either 'hier' or 'nodal', Bettina Schieche, 30th November 2011

d = size(nodelist,2); % dimension
basis = cell(d,1);

% find nodes in nodelist in all directions
if d == 1
    ind = (1:length(nodelist))';
    basis{1} = ind(levellist(ind) <= level);
else
	for i = 1:d
		I = [1:i-1,i+1:d]; % 1:d without i
		ind = find(max(abs(bsxfun(@minus,nodelist(:,I),node(:,I))),[],2) <= tol); % look for nodes in direction i
        if strcmp(flag,'hier')
            basis{i} = ind(levellist(ind,i) <= level(i)); % take only the nodes with level smaller or equal to the current level (hierarchical structure!)
        elseif strcmp(flag,'nodal')
            basis{i} = ind;
        end
	end
end
	
	
	
	
