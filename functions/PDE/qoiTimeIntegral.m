function q = qoiTimeIntegral(sol,mesh,freeNodes,a,b)

% calculate the time integral of a given solution for all nodes given by freeNodes
%
% Bettina Schieche, 2011

    quadTol = 1e-4;
    q = quadv(@(t)subfun(sol,mesh,freeNodes,t,a,b), ... 
              0,sol.x(end),quadTol);

end

function u = subfun(sol,mesh,freeNodes,t,a,b)

  u = zeros(mesh.nPoints,length(t));
  u(freeNodes,:) = deval(sol,t);
  u = mesh.valueAt(a,b,u);
  
end
