function y = parabolic2D(nodes,problemfolder,compute,qoiHandles)

% function to perform a numerical simulation for each collocation point
%
% nodes: collocation points for which to run the simulation (one point for each row)
% problemfolder: where to store the solution
% compute: function handle to start the simulation
% qoiHandles: quantities of interest to evaluate as function handles in a cell array
%
% Bettina Schieche, 2012

n = size(nodes,1); % number of nodes
allhelpvec = 1:n; % allhelpvec(i) should indicate which outfile number the nodes have
nodelistPath = [problemfolder,'/solutions/nodelist'];

if exist([nodelistPath,'.mat'],'file')==2  %if a list of nodes in solutions folder exists
    load(nodelistPath,'nodelist'); % load the list (each row = one point, last column number of kardos-outfile)
    oldnodes = nodelist(:,1:end-1);
    index = max(nodelist(:,end))+1; % next free outfile number
    vec = 1:n;
    for i = 1:n
        [J,I] = tolismember(nodes(i,:),oldnodes(:,1:size(nodes,2)),10^-8); % check if the node has been computed before
        if J==1
            vec(i) = 0; % indicates that nodes(i,:) needs no MATLAB run
            allhelpvec(i) = nodelist(I,end); % determine where to find the kardos solution for the current node
        else
            allhelpvec(i) = index;
            index = index + 1; % next free outfile number
        end
    end
    I = find(vec); % find the already existing nodes, marked with 0 before
    helpvec = allhelpvec(I); % outfile-numbers for kardos runs
    newnodes = nodes(I,:); % nodes for kardos runs
    n = size(newnodes,1); % number of truly new nodes
else % no solutions, yet
    nodelist = [];
    helpvec = allhelpvec;
    newnodes = nodes;
end

if n>0
    disp('Calculation started...')
end

solpath = [problemfolder,'/solutions/'];

for i = 1:n
    sol = compute(newnodes(i,:));
    sol = sol{:};
    if size(sol.y,1)<200000
        save([solpath,'solution',int2str(helpvec(i))],'-struct','sol');
    else
        sol.y = sol.y(73,:);
        sol.idata.dif3d = sol.idata.dif3d(73,:,:);
        savesol = sol;
        save([solpath,'solution',int2str(helpvec(i))],'-struct','savesol');
    end
    fid = fopen([solpath,'outfile',int2str(helpvec(i))],'w');
    for j = 1:length(qoiHandles)
        % Computation of the qoi via function handle.
        fprintf(fid,'%1.10e\n',qoiHandles{j}(sol));
    end
    fclose(fid);
end

if n >= 1
    nodelist = [nodelist;[newnodes,helpvec']]; % new nodelist
    save([nodelistPath,'.mat'],'nodelist') % save new nodelist
end

y = [];

% read QoI's
for i = 1:length(allhelpvec)
    newy = textread([solpath,'outfile',int2str(allhelpvec(i))])';
    y = [y;[newy newy.^2]];
end
end

