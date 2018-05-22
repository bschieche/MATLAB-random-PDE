function q = qoiSpaceIntegral(sol,mesh,freeNodes)

% calculate the space integral of a given solution on a given mesh
%
% Bettina Schieche, 2011

    q = spaceIntegral(sol.x(end),sol,mesh,freeNodes);

end
