function [deltagrid,intInfo,rule] = fillDeltaGrid(a,b,slevel,rule,d)

% deltagrid = structure with all rules that occur; each one is created by deltaGrid
% intInfo = cell array with as many entries as rules; each one is created by getIntInfo

if ~iscell(rule)
    rule = {rule,ones(1,d)};
end
deltagrid = [];
intInfo = cell(length(rule)-1,1);
for i = 1:length(rule)-1
    if strcmp(rule{i},'GaussPattersonLegendre') % other treatment of the 1D node structure
        deltagrid.(char(rule{i})) = deltaGrid(a(i),b(i),slevel,rule{i},1); % get the information about the quadrature rules Delta^i
        intInfo{i} = getIntInfo(deltagrid.(char(rule{i})),a(i),b(i),1);
    elseif strcmp(rule{i},'GaussPattersonHermite') % special treatment of the 1D node structure
        deltagrid.(char(rule{i})) = deltaGridGaussPattersonHermite(); % get the information about the quadrature rules Delta^i
        intInfo{i} = getIntInfoGaussHermitePatterson(deltagrid.(char(rule{i})));
    elseif strcmp(rule{i},'fclencurtgauss') % special treatment of the 1D node structure
        deltagrid.(char(rule{i})) = deltaGrid(a(i),b(i),slevel,rule{1}); % get the information about the quadrature rules Delta^i
        intInfo{i} = getIntInfoFclencurtgauss(deltagrid.(char(rule{i})));
    else
        deltagrid.(char(rule{i})) = deltaGrid(a(i),b(i),slevel,rule{i}); % get the information about the quadrature rules Delta^i
        intInfo{i} = getIntInfo(deltagrid.(char(rule{i})),a(i),b(i));
    end
end
