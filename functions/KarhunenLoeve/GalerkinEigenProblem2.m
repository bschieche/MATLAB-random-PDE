function [V,lambda] = GalerkinEigenProblem2(L,coordinates,elements,type)

% approximate the eigenvalues of a random process with covariance function of the Bessel type
% of the second kind and of first order by a Galerkin projection 

% you can get the eigenvalues for sigma~=1 by multiplying them by sigma^2
%
% Bettina Schieche, 2010

sigma = 1; %standard deviation of the process

ind1 = [elements(:,1) elements(:,1) elements(:,1) ...
        elements(:,2) elements(:,2) elements(:,2) ...
        elements(:,3) elements(:,3) elements(:,3)];

ind2 = [elements(:,1) elements(:,2) elements(:,3) ...
        elements(:,1) elements(:,2) elements(:,3) ...
        elements(:,1) elements(:,2) elements(:,3)];
    
ne = size(ind1,1);

x1 = coordinates(elements(:,1),1);
y1 = coordinates(elements(:,1),2);
x2 = coordinates(elements(:,2),1);
y2 = coordinates(elements(:,2),2);
x3 = coordinates(elements(:,3),1);
y3 = coordinates(elements(:,3),2);

Delta = (x2-x1).*(y3-y1)-(y2-y1).*(x3-x1);
volT = 0.5*abs(Delta); % |T| = 0.5*abs(Delta)

% p_i(x,y) = 1+a_i(x-x_i)+b_i(y-y_i)
a1 = (y2-y3)./Delta;
a2 = (y3-y1)./Delta;
a3 = (y1-y2)./Delta;
b1 = (x3-x2)./Delta;
b2 = (x1-x3)./Delta;
b3 = (x2-x1)./Delta;
Delta = abs(Delta);
% compute mass matrix
m = 1/24*[2*Delta 1*Delta 1*Delta 1*Delta 2*Delta 1*Delta 1*Delta 1*Delta 2*Delta];
M = sparse(ind1,ind2,m,size(coordinates,1),size(coordinates,1));

% quadratrue points
c1x = (x2-x1)/6+(x3-x1)/6+x1;
c2x = 2*(x2-x1)/3+(x3-x1)/6+x1;
c3x = (x2-x1)/6+2*(x3-x1)/3+x1;
c1y = (y2-y1)/6+(y3-y1)/6+y1;
c2y = 2*(y2-y1)/3+(y3-y1)/6+y1;
c3y = (y2-y1)/6+2*(y3-y1)/3+y1;

% evaluate basis in all quadrature points and multiply with weights
c1 = [(1+a1.*(c1x-x1)+b1.*(c1y-y1)).*volT;(1+a2.*(c1x-x2)+b2.*(c1y-y2)).*volT;(1+a3.*(c1x-x3)+b3.*(c1y-y3)).*volT]/3;
c2 = [(1+a1.*(c2x-x1)+b1.*(c2y-y1)).*volT;(1+a2.*(c2x-x2)+b2.*(c2y-y2)).*volT;(1+a3.*(c2x-x3)+b3.*(c2y-y3)).*volT]/3;
c3 = [(1+a1.*(c3x-x1)+b1.*(c3y-y1)).*volT;(1+a2.*(c3x-x2)+b2.*(c3y-y2)).*volT;(1+a3.*(c3x-x3)+b3.*(c3y-y3)).*volT]/3;
% for 3 Point quadrature we need 3 of such c matrices

clear m a1 a2 a3 b1 b2 b3 x1 x2 x3 y1 y2 y3 Delta ind1 ind2 volT;

elements = elements(:);
tic
covInfo = cell(6,1);
for i = 1:6
	covInfo{i} = zeros(ne);
end
for i = 1:ne
	covInfo{1}(:,i) = besselcov(sqrt((c1x(i)-c1x).^2+(c1y(i)-c1y).^2),L,sigma);
	covInfo{2}(:,i) = besselcov(sqrt((c1x(i)-c2x).^2+(c1y(i)-c2y).^2),L,sigma);
	covInfo{3}(:,i) = besselcov(sqrt((c1x(i)-c3x).^2+(c1y(i)-c3y).^2),L,sigma);
	covInfo{4}(:,i) = besselcov(sqrt((c2x(i)-c2x).^2+(c2y(i)-c2y).^2),L,sigma);
	covInfo{5}(:,i) = besselcov(sqrt((c2x(i)-c3x).^2+(c2y(i)-c3y).^2),L,sigma);
	covInfo{6}(:,i) = besselcov(sqrt((c3x(i)-c3x).^2+(c3y(i)-c3y).^2),L,sigma);
end
toc
disp('Covariance information complete')

clear c1x c1y c2x c2y c3x c3y;

C = zeros(size(coordinates,1));
I = cell(size(coordinates,1),1);
J = cell(size(coordinates,1),1);
for k = 1:size(coordinates,1)
	I{k} = find(elements==k);
	J{k} = mod(I{k}-1,ne)+1;
	if isempty(I{k})
		continue;
	end
    for l = 1:k
        if isempty(I{l})
			continue;
        end
		newentry = 	c1(I{k})'*covInfo{1}(J{k},J{l})*c1(I{l})+...
					c1(I{k})'*covInfo{2}(J{k},J{l})*c2(I{l})+...
					c1(I{k})'*covInfo{3}(J{k},J{l})*c3(I{l})+...
					c2(I{k})'*covInfo{2}(J{k},J{l})*c1(I{l})+...
					c2(I{k})'*covInfo{4}(J{k},J{l})*c2(I{l})+...
					c2(I{k})'*covInfo{5}(J{k},J{l})*c3(I{l})+...
					c3(I{k})'*covInfo{3}(J{k},J{l})*c1(I{l})+...
					c3(I{k})'*covInfo{5}(J{k},J{l})*c2(I{l})+...
					c3(I{k})'*covInfo{6}(J{k},J{l})*c3(I{l});
        C(k,l) = newentry;
        C(l,k) = newentry;
    end
end
disp('Matrices for Eigenvalue problem complete')
clear covInfo c1 c2 c3;
%%
tic
switch type
	case 'sparse'
		[V,D] = eigs(@(x)(C*x),size(M,1),M,min(200,size(M,1)));
	case 'full'
		[V,D] = eig(C,full(M));
end
toc
lambda = diag(D);
[lambda,I] = sort(lambda,'descend');
V = V(:,I);

fprintf('\n3-Point rule, case = %s\n',type)
fprintf('Eigenvalue check: sum(lambda) = %1.10e\n',sum(lambda));
fprintf('Orthonormality check: norm(V^t*M*V)= %1.10e\n\n',norm(V(:,1:10)'*M*V(:,1:10)));
