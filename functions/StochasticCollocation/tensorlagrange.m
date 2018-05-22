function y = tensorlagrange(x,nodelist,node,basis)

% y = tensorlagrange(x,nodelist,node,basis)
%
% function to evaluate tensor products of lagrange polynomials in x
%
% nodelist = interpolation nodes
% node = determines which basis function is wanted
% basis = cell array, each entry = vector of indices for nodes in contributing to the basis function
%
% Bettina Schieche, 27th of July 2010

d = size(node,2); % dimension

y = ones(size(x,1),1);

for i = 1:d
	y = y.*lagrange(x(:,i),nodelist(basis{i},i),node(i)); % tensor product of langrange polynomials
end

