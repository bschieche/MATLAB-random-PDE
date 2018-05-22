function q = qoiSpaceTimeIntegral(sol,mesh,freeNodes)

% calculate the space time integral of a given solution on a given mesh
%
% Bettina Schieche, 2011

    quadTol = 1e-4;
    q = quadv(@(t)spaceIntegral(t,sol,mesh,freeNodes),0,sol.x(end),quadTol);

end
