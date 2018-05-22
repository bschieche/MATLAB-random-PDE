function u = qoiMaxValue(sol,mesh,freeNodes,a,b)

% evaluate a solution on a given mesh in points given by freeNodes
%
% Bettina Schieche, 2011

  u = zeros(mesh.nPoints,1);
  u(freeNodes,:) = sol.y(:,end);
  u = mesh.valueAt(a,b,u);

end
