function [newsurplus,sol,interpval] = getSurplus(nodes,allnodes,level,func,surplus,problemfolder,rule)

% surplus = getSurplus(nodes,allnodes,level,func,surplus)
%
% function to calculate hierachical surplus for new nodes
%
% nodes = new nodes
% allnodes = old nodes
% level = level indices for each node in allnodes and each dimension
% func = function to interpolate (string)
% surplus = hierarchical surpluses which belong to allnodes
% 
% newsurplus = hierarchical surpluses which belong to nodes
%
% Bettina Schieche, 27th of July 2010
% returns also the solution itself now: sol, 2nd of August 2010
% also returns higher moments now: 4th of August 2010
% also returns interpval now: 5th of August 2010
% gets "problemfolder" as input: a string determining the KARDOS folder: 11th of October 2010
% gets "rule" in order to know which dimension belongs to which rule: 3nd of February 2001

newpointno = size(nodes,1); % number of new nodes
pointno = size(allnodes,1); % number of old nodes
interpval = zeros(newpointno,pointno);
for l = 1:pointno % evaluate older basis functions one by one
    % get information about nodes in all directions with level smaller or equal to the current level
    basis = node2basis(allnodes,allnodes(l,:),level,level(l,:),10^-8,'hier');
    % evaluate l. basis function in nodes
    y = tensorlagrange(nodes,allnodes,allnodes(l,:),basis);
    interpval(:,l) = y;
end
% calculate new hierarchical surpluses
sol = feval(func,nodes,problemfolder);
for i = 1:length(surplus)
    newsurplus{i} = sol.^i-interpval*surplus{i}; %compute higher moments
end
