function out = spaceIntegral(t,U,mesh,freeNodes)

% calculate the space integral for a given point in time t for a given solution U
%
% Bettina Schieche, 2011

  u = zeros(mesh.nPoints,length(t));
  u(freeNodes,:) = deval(t,U);
  out = sum(1/3*bsxfun(@times,u(mesh.cells(:,1),:) ...
                             +u(mesh.cells(:,2),:) ...
                             +u(mesh.cells(:,3),:),mesh.cellIntegral));

