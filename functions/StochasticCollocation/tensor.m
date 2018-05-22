function [B,C] = tensor(A,W,L)
% function to find all combinations of nodes that result in a tensor product of 1D-quadrature rules, it is called
% recursively!
%
% A: cell array containing vectors of nodes
% W: cell array containing vectors of weights (same length as the vectors in A!)
% L: vector of the same length as A and W, determining for each vector in A (resp. in W) how long it is, namely how many
% nodes (resp. weights) there are in it
%
% B = matrix containing all possible combinations of nodes (one row = one node in the space of multiple dimensions)
% C = matrix containing the corresponding weights (all entries in one row have to be multiplied after running this
% function in order to get the weights for the quadrature rule in multiple dimensions)
%
% Bettina Schieche, 11th of September 2009

A1 = A{1}; % get the first vector of nodes (can be seen as the nodes for the first coordinate in multiple dimensions)
A2 = {A{2:end}}; % get A without the first entry
W1 = W{1}; % same like A1 for weights now
W2 = {W{2:end}}; % same like A2 for weights now
L1 = L(2:end); % get L without the length corresponding to the first vetor of nodes

if isempty(A2) == 1
    % if A contains just one vector just A1 and W1 are returned
    B = A1;
    C = W1;
    return
else
    % if there is more than one vector of nodes, A1 is replicated and sorted 
    [B,I] = sort(repmat(A1,prod(L1),1));
    % the weights are rearranged in the same order like B
    C = repmat(W1,prod(L1),1);
    C = C(I);
    [B1,C1] = tensor(A2,W2,L1); % call tensor by forgetting the information corresponding to first dimension
    B = [B,repmat(B1,length(A1),1)]; % replicate the output of tensor and add it to B and C
    C = [C,repmat(C1,length(A1),1)];
end
        


        