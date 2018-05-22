function deltagrid = deltaGrid(a,b,k,rule,varargin)
% function that returns all information needed for the 1D-quadrature rules Delta^i for the special case of
% Clenshaw-Curtis rule (see papers for notation explanation)
%
% [a,b]: 1D-interval, considered range
% d: dimension, needed in order to decide how many rules there have to be constructed
% k: Smolyak level
%
% deltagrid: cell array, first column: nodes, second column: weights, third column: no. of nodes resp. weights
% all nodes and weights are column vectors
%
% Bettina Schieche, 11th September 2009
% changed to rule-depending function, Bettina Schieche, 22th January 2010

m = 2.^(0:k+1)+1; % calculate order of the quadrature rules (necessary for nested grids!)
m(1) = 1; % set the first one manually

% should be called if the endpoints are not nodes themselves (open rules). In this case the nestedness behaves slightly different
if nargin > 4
    m = 2.^(1:k+2)-1;
end
deltagrid = cell(k+2,3); % initialize deltagrid

[deltagrid{1,1},deltagrid{1,2}] = feval(rule,m(1),a,b);  % compute the first order quadrature rule
deltagrid{1,3} = 1; % the first order rule contains of one node

wbefore = deltagrid{1,2}; % needed for later comparison
for j = 2:length(m)
    [x,w] = feval(rule,m(j),a,b); % get the nodes and weights for order m(j)
    deltagrid{j,1} = x; % save the nodes
    w1 = w; % keep the nodes for later usage
    % find the indices of x, whose nodes appear in xbefore as well (j = 2 is a special case!)
    if nargin > 4
        I = 2:2:length(x);
    else
        if j==2
            I = 2;
        else
            I = 1:2:length(x);
        end
    end
    
    % compute the weights, Delta^j consists of the same nodes like x, but the weights corresponding to xbefore have to
    % be subtracted
    w1(I) = w1(I)-wbefore;
    deltagrid{j,2} = w1;
    deltagrid{j,3} = length(x);
    % save the current nodes and weights in order to compare them again with x and w computed in the next loop
    wbefore = w;
end





