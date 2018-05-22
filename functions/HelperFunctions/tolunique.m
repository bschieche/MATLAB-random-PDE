function [y,I,J] = tolunique(x,tol)

% function [y,I,J] = tolunique(x,tol)
% 
% find the unique elements of x with a tolerance tol
% x must be sorted before!!!!!!!!!!
% 
% y: unique elements
% I: y = x(I)
% J: x = y(J)
%
% Bettina Schieche, 19th of March 2010

n = size(x,1); % number of elements to sort
z = ones(n,1); % vector for checking purposes

j = 1;
k = 1;
J = zeros(n,1);
J(1) = 1;
% step through x: if the distance between 2 elements is smaller than tol, it is marked by zero in z
for i = 2:n
    if sum(abs(x(i,:)-x(j,:)))<=tol
        z(i) = 0;
        J(i) = k;
    else
        j = i;
        k = k+1;
        J(i) = k;
    end
end
I = find(z);
y = x(I,:); % extract the unique elements 
