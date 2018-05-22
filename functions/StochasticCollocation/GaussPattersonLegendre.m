function [x,w] = GaussPattersonLegendre(N1,a,b)

% The following is only exact for 7-8 digits and I don't understand the
% exact meaning of N1
%[x,w] = nwspgr('KPU',1,N1);
%x = (x-0.5)*(b-a)+(b+a)/2;
%w = w*(b-a);

[x,w] = patterson_rule(N1,a,b);
