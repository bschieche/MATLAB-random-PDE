function y = besselcov(r,b,sigma)
% evaluates a bessel covariance function with correlation length b and
% standard deviation sigma
%
% Bettina Schieche, 2010

y=r/b.*besselk(1,r/b)*sigma^2;
I=find(r==0);
y(I)=sigma^2;
