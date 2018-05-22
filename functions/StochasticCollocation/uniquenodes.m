function [nodes,varargout] = uniquenodes(nodes,varargin)

% function to get the unique nodes and sum up the corresponding weights
%
% nodes: d columns, each row is one node
% weights: vector the same length as nodes
%
% Bettina Schieche, 22th January 2010
%
% can also sort the level information, Bettina Schieche, 30th November 2011

[nodes,I] = sortrows(nodes); % sort the nodes
[nodes,K,J] = tolunique(nodes,10^-8); % get the unique nodes

if nargin >= 2 % contains one vector of weights
    weights = varargin{1}(I); % bring the weights in the same order
    uniqueWeights = zeros(size(nodes,1),1); % initialize the new weights vector
    for j = 1:max(J)
        uniqueWeights(j) = sum(weights(J==j)); % sum up the values that belong to one node
    end
    varargout{1} = uniqueWeights;
end

if nargin >= 3 % contains level matrix
    level = varargin{2}(I,:); % bring the level in the same order
    varargout{2} = level(K,:);
end