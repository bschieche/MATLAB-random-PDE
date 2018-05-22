function ind = findActiveNodes(active,level)

% helper function to find the currently active interpolation nodes
%
% Bettina Schieche, 2010

ind = [];
for i = 1:size(active,1)
    ind = [ind;find(sum(abs(bsxfun(@minus,level,active(i,:))),2)==0)];
end
