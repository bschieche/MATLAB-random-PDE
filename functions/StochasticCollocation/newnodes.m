function [nodes,varargout] = newnodes(deltagrid,j,allnodes,rule)

% nodes = newnodes(deltagrid,j,allnodes,rule)
%
% function to add nodes according to index j
%
% deltagrid = information about nodes on different levels, see deltaGrid.m for further information
% j = index in Smolyak sum
% allnodes = oldnodes
%
% nodes = new nodes
%
% Bettina Schieche, 27th of July 2010
%
% updated deltagrid, which is now a structure with more than one deltagrids
% deltagrid is now a structure of delatgrids, rule is needed to know which
% dimension belongs to which rule, 2nd of February 2011
%
% if allnodes is empty, a whole new Smolyak rule is calculated based on
%   some level stored in j; the weights are then stored in varargout, 2nd of
%   November 2011
%
% is now able to return the corresponding level information, 30th November 2011

nodes = [];
weights = [];
level = [];
for i = 1:size(j,1)
    [tensorinfo_vec,tensorinfo_w,tensorinfo_length] = tensorinfo(rule,deltagrid,j(i,:));
    [newNodes,newWeights] = tensor(tensorinfo_vec,tensorinfo_w,tensorinfo_length);
    nodes = [nodes;newNodes]; % get the nodes belonging to index j
    weights = [weights;newWeights];
    level = [level;repmat(j(i,:),size(newNodes,1),1)];
end
if isempty(nodes)
    weights = [];
    level = [];
else
    [nodes,weights,level] = uniquenodes(nodes,prod(weights,2),level); % get the unique nodes
end

if isempty(allnodes)
    varargout{1} = weights;
    varargout{2} = level;
    return
else
    newpointno = size(nodes,1); % number of new nodes
    check = zeros(newpointno,1); % to mark the truly new nodes
    for i = 1:newpointno
        [c,~] = tolismember(nodes(i,:),allnodes,10^-8); % check if the node has been computed before
        if ~c % if it is a truly new node
            check(i) = 1;
        end
    end
    nodes = nodes(check == 1,:); % reduce to truly new nodes
end