function [points,triangles,boundaryEdges,M1,M2] = RefineTriangularGrid(points,triangles,boundaryEdges)

  nT = size(triangles,1);
  nP = size(points,1);

  edges = [triangles(:,1) triangles(:,2);
           triangles(:,2) triangles(:,3);
           triangles(:,3) triangles(:,1)];
  toSwap = edges(:,1)>edges(:,2);
  edges(toSwap,[1 2]) = edges(toSwap,[2 1]);
  edges = unique(edges,'rows');

  nE = size(edges,1);

  % The entry a_ij of the matrix is 1, if the point i belongs to
  % edge j.
  i = edges(:)';
  j = [1:nE 1:nE]';
  s = ones(size(i));
  pointEdgeMatrix = sparse(i,j,s,nP,nE);

  % The entry a_ij of the matrix is 1, if the point i belongs to
  % triangle j.
  i = triangles(:)';
  j = [1:nT 1:nT 1:nT]';
  s = ones(size(i));
  pointTriangleMatrix = sparse(i,j,s,nP,nT);

  % The entry a_ij of the matrix is 1, if the edge i belongs to
  % triangle j.
  edgeTriangleMatrix = (pointEdgeMatrix' * pointTriangleMatrix) == 2;
  
  newPoints = [0.5*(points(edges(:,1),1)+points(edges(:,2),1)) ...
               0.5*(points(edges(:,1),2)+points(edges(:,2),2))];
             
  newPointsIndex = nP+(1:nE)';
  
  [edgesOfTriangle,~] = find(edgeTriangleMatrix);
  edgesOfTriangle = reshape(edgesOfTriangle,3,nT)';
  
  newTriangles = cell(4,1);
  
  for i=1:3
  
    pick23 = edges(edgesOfTriangle(:,1),1) ~= triangles(:,i) ...
           & edges(edgesOfTriangle(:,1),2) ~= triangles(:,i);
    pick31 = edges(edgesOfTriangle(:,2),1) ~= triangles(:,i) ...
           & edges(edgesOfTriangle(:,2),2) ~= triangles(:,i);
    pick12 = edges(edgesOfTriangle(:,3),1) ~= triangles(:,i) ...
           & edges(edgesOfTriangle(:,3),2) ~= triangles(:,i);

    assert(isequal(pick23,~(pick31|pick12)));
    assert(isequal(pick31,~(pick12|pick23)));
    assert(isequal(pick12,~(pick23|pick31)));

    newTriangles{i} = nan(nT,3);
    
    newTriangles{i}(pick23,:) = [                     triangles(pick23,i) ...
                              newPointsIndex(edgesOfTriangle(pick23,2)) ...
                              newPointsIndex(edgesOfTriangle(pick23,3))];
    newTriangles{i}(pick31,:) = [                     triangles(pick31,i) ...
                              newPointsIndex(edgesOfTriangle(pick31,3)) ...
                              newPointsIndex(edgesOfTriangle(pick31,1))];
    newTriangles{i}(pick12,:) = [                     triangles(pick12,i) ...
                              newPointsIndex(edgesOfTriangle(pick12,1)) ...
                              newPointsIndex(edgesOfTriangle(pick12,2))];
                             
  end
  
  newTriangles{4} = [newPointsIndex(edgesOfTriangle(:,1)) ...
                     newPointsIndex(edgesOfTriangle(:,2)) ...
                     newPointsIndex(edgesOfTriangle(:,3))];
  
  points = [points; newPoints];
  
  triangles = nan(4*size(triangles,1),3);
  triangles(1:4:end,:) = newTriangles{1};
  triangles(2:4:end,:) = newTriangles{2};
  triangles(3:4:end,:) = newTriangles{3};
  triangles(4:4:end,:) = newTriangles{4};
  
  if nargin >= 3 && nargout >= 3
  
    newBoundaryEdges = cell(2,1);  
    globalEdgeNumber = pointEdgeMatrix(boundaryEdges(:,1),:) ...
                     & pointEdgeMatrix(boundaryEdges(:,2),:);
    [globalEdgeNumber,~] = find(globalEdgeNumber');
    newBoundaryEdges{1} = [boundaryEdges(:,1) newPointsIndex(globalEdgeNumber)];
    newBoundaryEdges{2} = [newPointsIndex(globalEdgeNumber) boundaryEdges(:,2)];

    boundaryEdges = nan(2*size(boundaryEdges,1),2);
    boundaryEdges(1:2:end,:) = newBoundaryEdges{1};
    boundaryEdges(2:2:end,:) = newBoundaryEdges{2};
  
  else
    
    boundaryEdges = [];
    
  end
  
  if nargout >=4
    
    i = [(1:nP)'; newPointsIndex; newPointsIndex];
    j = [(1:nP)'; edges(:,1); edges(:,2)];
    s = [ones(nP,1); 0.5*ones(2*nE,1)];
    M1 = sparse(i,j,s,size(points,1),nP);
    
  else
    
    M1 = [];
    
  end
  
  if nargout >=5
    
    i = (1:nP)';
    j = (1:nP)';
    s = ones(nP,1);
    M2 = sparse(i,j,s,nP,size(points,1));
    
  else
    
    M2 = [];
    
  end
  
  
end