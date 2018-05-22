function [allnodes,level,surplus,intval,nodesol,alpha,deltagrid] = ...
    anisoHierLagrangeSimplified(a,b,slevel,rule,d,func,tol,omega,problemfolder,maxcounter,qoi)

% function to interpolate a function 'func' in d dimensions on [a,b]^d
% using maximal level 'slevel'. The generalized Smolyak algorithm is used
% with global hierarchical Lagrange polynomials as basis functions
%
% rule = string to determine the quadrature nodes which form the grid
% func = string to determine the underlying functions which needs to be evaluated in the nodes
% tol = tolerance for the global error indicator
% omega = balancing number between 0 and 1 (the higher the stronger anisotropic)
%
% allnodes = resulting nodes (each row = one node)
% level = the same size as allnodes with level indices for each node and each dimension
% surplus = weight for each node
%
% Bettina Schieche, 21th of July 2010
% surplus changed to a cell array of hierarchical surpluses for moments of the function -> input parameter 'moments'
% also returns nodesol and weights when resorting the sum to a linear combination of function evaluations and corresponding weights
% gets "name" as input: a string determining the folder for solutions, 14th of September 2010
% gets "problemfolder" as input: a string determining the main solution folder, 11th of October 2010
% rule is now a cell array of several rules, where an additional entry is a vector which specifies which rule belongs to which dimension
% (e.g. [1 1 2 2] for d = 4 and 2 rules), 2nd of February 2011

% get the information about all difference 1D quadrature rules and all possible integrals of 1D-Lagrange polynomials
[deltagrid,intInfo,rule] = fillDeltaGrid(a,b,slevel,rule,d);
Id = eye(d); % identity matrix where unit vectors are taken from below
ind = ones(1,d); % starting index
inactive = []; % initialize inactive index set with the empty set
active = ind; % initialize the active index set with i

% get input for tensor.m
[tensorinfo_vec,tensorinfo_w,tensorinfo_length] = tensorinfo(rule,deltagrid,ind);
% compute the nodes corresponding to the first index
allnodes = tensor(tensorinfo_vec,tensorinfo_w,tensorinfo_length);
level = ind; % initialize level
nodes = newnodes(deltagrid,repmat(ind,d,1)+Id,allnodes,rule);
feval(func,[allnodes;nodes],problemfolder);
surplus{1} = feval(func,allnodes,problemfolder); % evaluate first node
surplus{2} = surplus{1}.^2; % get the square of the first node
nodesol = surplus{1}; % save solution of first node
intval = 1; % integral across first basis functions
alpha = intval; % weight corresponding to the first node

g = 1; % first local error indicator
eta = g;
allErr = cell(d,1);
counter = 0;

while eta > tol &&  ~isempty(g)
    counter = counter + 1;
    if counter > maxcounter
        break
    end
    [active,inactive,g,~,ind] = updateInd(active,inactive,g,eta,omega);
    collectj = [];
    dimind = zeros(d,1);
    for k = 1:d
        j = ind + Id(k,:); % add the k. unit vector
        R = repmat(j,d,1); % replicate j in order to check whether it is an admissible index (see the rows below)
        % find the rows where the diagonal element is not 1 (because the admissibility condition has to be checked for elements less than 1)
        I = find(diag(R) > 1);
        % check if j is admissible
        if all(ismember(R(I,:)-Id(I,:),inactive,'rows'))
            collectj = [collectj;j];
            dimind(k) = 1;
        end
    end
    nodes = newnodes(deltagrid,collectj,allnodes,rule); % maybe nodes has another order than the nodes which are saved in allnodes later (just as a note!)
    if isempty(nodes)
        errorInfo(0,max(g),eta);
        continue
    end
    feval(func,nodes,problemfolder);
    for k = find(dimind == 1)'
        j = ind + Id(k,:); % add the k. unit vector
        R = repmat(j,d,1); % replicate j in order to check whether it is an admissible index (see the rows below)
        % find the rows where the diagonal element is not 1 (because the admissibility condition has to be
        % checked for elements less than 1)
        I = find(diag(R) > 1);
        % check if j is admissible
        if all(ismember(R(I,:)-Id(I,:),inactive,'rows'))
            nodes = newnodes(deltagrid,j,allnodes,rule);
            if isempty(nodes)
                continue
            end
            active = [active;j]; % add j to the active index set
            newintval = getIntval(nodes,intInfo,d,rule);
            intval = [intval;newintval];
            [newsurplus,newnodesol,interpval] = getSurplus(nodes,allnodes,level,func,surplus,problemfolder,rule);
            for i = 1:length(surplus)
                surplus{i} = [surplus{i};newsurplus{i}];
            end
            nodesol = [nodesol;newnodesol];
            alpha = [[alpha zeros(size(allnodes,1),size(nodes,1))];[-interpval*alpha,eye(size(nodes,1))]];
            level = [level;repmat(j,size(nodes,1),1)];
            allnodes = [allnodes;nodes];
            allErr{k} = [allErr{k};localError(newsurplus,newintval,surplus,intval,qoi)]; % add local error indicator
            if abs(newintval'*newsurplus{1}(:,qoi)) < 0
                active = active(1:end-1,:);
            else
                g = [g;localError(newsurplus,newintval,surplus,intval,qoi)]; % add local error indicator
            end
        end
    end
    I = findActiveNodes(active,level);
    % type of error indicator: standard deviation
    weights = intval'*alpha;
    integmean = weights*nodesol(:,qoi);
    integstd = sqrt(weights*nodesol(:,qoi).^2-integmean.^2);
    J = setdiff(1:length(weights),I);
    weights2 = intval(J)'*alpha(J,J);
    integmean2 = weights2*nodesol(J,qoi);
    integstd2 = sqrt(weights2*nodesol(J,qoi).^2-integmean2.^2);
    eta = max(abs((integstd-integstd2)./integstd));
    
    errorInfo(0,max(abs(g)),eta);
end
disp('end of loop');
