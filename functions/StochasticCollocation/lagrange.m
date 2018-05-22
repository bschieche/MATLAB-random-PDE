function y = lagrange(x,nodelist,node)

% y = lagrange(x,nodelist,node)
%
% function to evaluate a 1D lagrange polynomial in x
%
% nodelist = interpolation nodes
% node = determines lagrange polynomial
%
% Bettina Schieche, 27th of July 2010

I = find(abs(nodelist-node)<=10^-10); % find node in nodelist
if isempty(I) || length(I) > 1 % throw error if node is not or multiple times in nodelist
	disp('error occurred');
	return;
end

y = ones(size(x));
for i = 1:size(nodelist,1)
	if i~=I
		y = y.*(x-nodelist(i))./(nodelist(I)-nodelist(i));
	end
end
