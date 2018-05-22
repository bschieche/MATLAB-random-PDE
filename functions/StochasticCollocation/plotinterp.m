function interpval = plotinterp(a,b,allnodes,level,surplus,nodes,flag,varargin)

% plotinterp(a,b,allnodes,level,surplus)
%
% procedure to plot the interpolant, only possible for 2 dimensions
%
% [a,b]^2 = domain
% allnodes = interpolation nodes
% level = level information for all nodes
% surplus = hierarchical surpluses for all nodes
%
% Bettina Schieche, 27th of July 2010
%
% Bettina Schieche, 1st of September 2011: varargin for a dimchecker: vector of 2 components specifying which
% dimensions to plot, the rest is filled with zeros, example:  plotinterp([-7 -7],[7 7],allnodes,level,surplus{1}(:,1),[],[7 8]);
%
% Bettina Schieche, 30th of November 2011: is now able to interpolate several surplus column vectors in one matrix simultaneously
% and gets a flag reading either 'hier' or 'nodal'
%
% Bettina Schieche, 7th of March 2012: for pure interpolation (no plotting) varargin{1} can be an index; 
% then the interpolation starts with this index instead of 1 (see lstart below)

if isempty(nodes)
    x = a(1):(b(1)-a(1))/100:b(1);
    y = a(2):(b(2)-a(2))/100:b(2);
    [X,Y] = meshgrid(x,y);
    lstart = 1;
else
    if nargin == 8 && ~isempty(varargin{1})
        lstart = varargin{1};
    else
        lstart = 1;
    end
end

interpval = 0;
for l = lstart:size(allnodes,1);
    basis = node2basis(allnodes,allnodes(l,:),level,level(l,:),10^-8,flag);
    if isempty(nodes)
        if isempty(varargin)
            interpval = interpval + surplus(l)*...
                tensorlagrange([reshape(X,length(x)^2,1) reshape(Y,length(x)^2,1)],allnodes,allnodes(l,:),basis);
        else
            dimchecker = varargin{1};
            vec = zeros(length(x)^2,size(allnodes,2));
            vec(:,dimchecker(1)) = reshape(X,length(x)^2,1);
            vec(:,dimchecker(2)) = reshape(Y,length(x)^2,1);
            interpval = interpval + surplus(l)*tensorlagrange(vec,allnodes,allnodes(l,:),basis);
        end
    else
        interpval = interpval + bsxfun(@times,surplus(l,:),tensorlagrange(nodes,allnodes,allnodes(l,:),basis));
    end
end

if isempty(nodes)
    mesh(X,Y,reshape(interpval,length(x),length(x)),'FaceColor','interp','EdgeColor','interp');
end
