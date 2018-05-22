function [active,inactive,g,eta,ind] = updateInd(active,inactive,g,eta,omega)

% [active,inactive,g,eta,ind] = updateInd(active,inactive,g,eta,omega)
%
% function to choose next index and update nodes and error indicators
%
% active = active indices
% inactive = inactive indices
% g = vector of local error indicators, the same length as active
% eta = global error indicator
% omega = balancing number between 0 and 1 (the higher the more anisotropic)
%
% Bettina Schieche, 27th of July 2010

[~,I] = min(sum(active,2)); %chose the point, whose level is the smallest (isotropic behaviour)
I = I(1); % w.l.o.g. take the first index
maxg = g(I); % get the corresponding local error indicator
disp(num2str(sum(active(I,:))));
disp(num2str((1-omega)*max(sum([active;inactive],2))));
if sum(active(I,:)) > (1-omega)*max(sum([active;inactive],2))
    [maxg,I] = max(abs(g)); % find maximal error indicator and the corresponding indices (anisotropic behaviour)
    I = I(1); % take w.l.o.g. the first found index
	disp('Anisotropic behaviour!');
else
	disp('Isotropic behaviour!');
end
J = setdiff(1:length(g),I); % find the complement
ind = active(I,:); % get the corresponding index
active = active(J,:); % delete the index from the active set
g = g(J); % delete the corresponding local error indicator
inactive = [inactive;ind]; % add index to the inactive index set
eta = eta-maxg; % subtract maxg from the global error estimator
