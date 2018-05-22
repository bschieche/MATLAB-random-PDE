classdef LinearTriangularMesh < LinearMesh
% Store triangular mesh data and perform mesh-related computations.
%
% Sebastian Ullmann
% Numerical Analysis and Scientific Computing
% Technische Universität Darmstadt
% 29/12/2010
% 
% See also LinearTriangularMesh/LinearTriangularMesh, LinearMesh.

  methods (Access = public)
    
    function obj = LinearTriangularMesh(points,triangles)
    %LinearTriangularMesh Create an instance of LinearTriangularMesh.
    %
    %Syntax:
    %
    %  obj = LinearTriangularMesh(points,tetrahedra)
    %
    %Arguments:
    %
    %  points       An array of mesh points.
    %
    %  triangles    An array of triangles. (Indices into points.)
    %
    %  obj          The LinearTriangularMesh object.
    %
    %See also LinearMesh.
    %
      nDim = 2;
      obj  = obj@LinearMesh(points,triangles,nDim);
    end
    function out = valueAt(obj,xp,yp,u)
    %valueAt Find value at point by interpolation.
        
      assert(all(size(xp)==size(yp)), ...
             'Arguments xp and yp must be of equal size.');
      assert(length(xp) == numel(xp), ...
             'Arguments xp and yp must be scalars or vectors.');
      assert(size(u,1) == obj.nPoints);
           
      A = obj.points(obj.cells(:,1),:);
      B = obj.points(obj.cells(:,2),:);
      C = obj.points(obj.cells(:,3),:);

      DeltaP = (B(:,1)-A(:,1)).*(C(:,2)-A(:,2)) ...
             - (B(:,2)-A(:,2)).*(C(:,1)-A(:,1));

      if iscell(u)
        out = cell(size(u));
        for i=1:numel(u)
          out{i} = zeros(length(xp),size(u{1},2));
        end
      elseif isnumeric(u)
        out = zeros(length(xp),size(u,2));
      else
        error('input u must be cell or numeric');
      end
        
      for i=1:length(xp)
        
        P1 = xp(i);
        P2 = yp(i);

        DeltaA = (B(:,1)-P1).*(C(:,2)-P2) ...
               - (B(:,2)-P2).*(C(:,1)-P1);

        DeltaB = (C(:,1)-P1).*(A(:,2)-P2) ...
               - (C(:,2)-P2).*(A(:,1)-P1);

        DeltaC = (A(:,1)-P1).*(B(:,2)-P2) ...
               - (A(:,2)-P2).*(B(:,1)-P1);

        DeltaMat = [DeltaA./DeltaP ...
                    DeltaB./DeltaP ...
                    DeltaC./DeltaP];

        idx = find(sum(DeltaMat>=0,2)==3,1);
        
        if iscell(u)
          for j=1:numel(u)
            out{j}(i,:) = DeltaA(idx)/DeltaP(idx).*u{j}(obj.cells(idx,1),:) ...
                        + DeltaB(idx)/DeltaP(idx).*u{j}(obj.cells(idx,2),:) ...
                        + DeltaC(idx)/DeltaP(idx).*u{j}(obj.cells(idx,3),:);
          end
        else
          out(i,:) = DeltaA(idx)/DeltaP(idx).*u(obj.cells(idx,1),:) ...
                   + DeltaB(idx)/DeltaP(idx).*u(obj.cells(idx,2),:) ...
                   + DeltaC(idx)/DeltaP(idx).*u(obj.cells(idx,3),:);
        end
        
      end
        
    end  
    function [bEdges bNormal bTrianglesIndex] = boundary(obj)
      
      nT = size(obj.cells,1);
      nP = size(obj.points,1);
    
      edges = [obj.cells(:,1) obj.cells(:,2);
               obj.cells(:,2) obj.cells(:,3);
               obj.cells(:,3) obj.cells(:,1)];
      toSwap = edges(:,1)>edges(:,2);
      edges(toSwap,[1 2]) = edges(toSwap,[2 1]);
      edges = unique(edges,'rows');
      
      nE = size(edges,1);

      % The entry a_ij of the matrix is 1, if the point i belongs to
      % edge j.
      i = edges(:)';
      j = [1:nE 1:nE]';
      s = ones(size(i));
      edgeMatrix = sparse(i,j,s,nP,nE);
      
      % The entry a_ij of the matrix is 1, if the point i belongs to
      % triangle j.
      i = obj.cells(:)';
      j = [1:nT 1:nT 1:nT]';
      s = ones(size(i));
      triangleMatrix = sparse(i,j,s,nP,nT);
     
      % The entry a_ij of the matrix is 1, if the edge i belongs to
      % triangle j.
      edgeTriangleMatrix = (edgeMatrix' * triangleMatrix) == 2;
      
      % Edges that are only part of one triangle. 
      bEdgesIndex = find(sum(edgeTriangleMatrix,2)==1);
      bEdges      = edges(bEdgesIndex,:);
      
      % Triangles that belong to the boundary Edges.
      [bTrianglesIndex,dummy] = find(edgeTriangleMatrix(bEdgesIndex,:)');
      bTriangles          = obj.cells(bTrianglesIndex,:);

      pick1 = bTriangles(:,1) ~= bEdges(:,1) ...
            & bTriangles(:,1) ~= bEdges(:,2);       
      pick2 = bTriangles(:,2) ~= bEdges(:,1) ...
            & bTriangles(:,2) ~= bEdges(:,2);
      pick3 = bTriangles(:,3) ~= bEdges(:,1) ...
            & bTriangles(:,3) ~= bEdges(:,2);
      
      thirdPoint = zeros(length(bEdgesIndex),1);
      thirdPoint(pick1) = bTriangles(pick1,1);
      thirdPoint(pick2) = bTriangles(pick2,2);
      thirdPoint(pick3) = bTriangles(pick3,3);
            
      a = [obj.points(bEdges(:,1),1)-obj.points(thirdPoint,1) ...
           obj.points(bEdges(:,1),2)-obj.points(thirdPoint,2) ...
           zeros(size(bEdges,1),1)];
               
      b = [obj.points(bEdges(:,2),1)-obj.points(thirdPoint,1) ...
           obj.points(bEdges(:,2),2)-obj.points(thirdPoint,2) ...
           zeros(size(bEdges,1),1)];
         
      e = [obj.points(bEdges(:,2),1)-obj.points(bEdges(:,1),1) ...
           obj.points(bEdges(:,2),2)-obj.points(bEdges(:,1),2) ...
           zeros(size(bEdges,1),1)];
      
      % Make bNormal perpendicular to the edges and with right length.
      bNormal = cross(a,b,2);
      bNormal = bsxfun(@rdivide,bNormal,sqrt(dot(bNormal,bNormal,2)));
      bNormal = cross(e,bNormal,2);
      % Make bNormal point outward.
      bNormal = bsxfun(@times,sign(dot(bNormal,a+b,2)),bNormal);
      bNormal = bNormal(:,[1 2]);
      
    end 
    function M = interpolationMatrix(obj,mesh)
      
      t1 = mesh.cells(:,1);
      t2 = mesh.cells(:,2);
      t3 = mesh.cells(:,3);
      
      A1 = mesh.points(t1,1);
      B1 = mesh.points(t2,1);
      C1 = mesh.points(t3,1);
      A2 = mesh.points(t1,2);
      B2 = mesh.points(t2,2);
      C2 = mesh.points(t3,2);
      
      ABAC = (B1-A1).*(C2-A2) - (B2-A2).*(C1-A1);
      signABAC = sign(ABAC);
      
      i = zeros(size(obj.points,1),3);
      j = zeros(size(obj.points,1),3);
      s = zeros(size(obj.points,1),3);

      for ip=1:size(obj.points,1)
        
        P1 = obj.points(ip,1);
        P2 = obj.points(ip,2);
        
        PB1 = B1-P1;
        PA1 = A1-P1;
        PC1 = C1-P1;
        PB2 = B2-P2;
        PA2 = A2-P2;
        PC2 = C2-P2;
       
        idx = find((PB1.*PC2 - PB2.*PC1).*signABAC >= 0 ...
                 & (PC1.*PA2 - PC2.*PA1).*signABAC >= 0 ...
                 & (PA1.*PB2 - PA2.*PB1).*signABAC >= 0, 1, 'first');
        
        PBPC = PB1(idx).*PC2(idx) - PB2(idx).*PC1(idx);
        PCPA = PC1(idx).*PA2(idx) - PC2(idx).*PA1(idx);          
        PAPB = PA1(idx).*PB2(idx) - PA2(idx).*PB1(idx);
        
        i(ip,:) = ip;
        j(ip,:) = [t1(idx) t2(idx) t3(idx)];
        s(ip,1) = PBPC/ABAC(idx);
        s(ip,2) = PCPA/ABAC(idx);
        s(ip,3) = PAPB/ABAC(idx);
       
      end    
      
      M = sparse(i,j,s,size(obj.points,1),size(mesh.points,1));
      
    end
    
  end
  
  methods (Access = protected)
    
    function out = computeDdx(obj,dir)
                
      if isempty(who('dir'))
        out = {obj.computeDdx(1) obj.computeDdx(2)};
        return
      end
      
      assert(isnumeric(dir) && isscalar(dir) && any(dir==[1 2]), ...
             'Choose 1 or 2 as input parameter.');
           
      Delta = obj.mixedProduct;

      % The vector ddxi is the x-derivative of the hat function that is
      % equal to 1 in the ith corner of the tetrahedron. Equivalently, if
      % phi_i is the hat function that is equal to 1 in the ith corner of
      % the tetrahedron and phi_i = a_i + b_i*x + c_i*y + d_i*z within the
      % tetrahedron, then ddxi = b_i, ddyi = c_i and ddzi = d_i.
          
      switch dir
        
        case 1
          
          ddx1 = (obj.points(obj.cells(:,2),2)-obj.points(obj.cells(:,3),2))./Delta;
          ddx2 = (obj.points(obj.cells(:,3),2)-obj.points(obj.cells(:,1),2))./Delta;
          ddx3 = -ddx1 -ddx2;
          
          
        case 2
          
          ddx1 = (obj.points(obj.cells(:,3),1)-obj.points(obj.cells(:,2),1))./Delta;
          ddx2 = (obj.points(obj.cells(:,1),1)-obj.points(obj.cells(:,3),1))./Delta;
          ddx3 = -ddx1 -ddx2;
          
        otherwise
          
          error('Choose dir as 1 or 2.');
          
      end
                
      % This is where to put the entries of the element matrix.
      i = [1:size(obj.cells,1) 1:size(obj.cells,1) 1:size(obj.cells,1)]';
      j = [obj.cells(:,1); obj.cells(:,2); obj.cells(:,3)];

      % Matrix assembly.
      out = sparse(i, j, [ddx1;ddx2;ddx3], size(obj.cells,1), size(obj.points,1), length(i));
    end
    function out = computeMixedProduct(obj)
      out = (obj.points(obj.cells(:,2),1)-obj.points(obj.cells(:,1),1)) ...
          .*(obj.points(obj.cells(:,3),2)-obj.points(obj.cells(:,2),2)) ...
          - (obj.points(obj.cells(:,2),2)-obj.points(obj.cells(:,1),2)) ...
          .*(obj.points(obj.cells(:,3),1)-obj.points(obj.cells(:,2),1));
    end
    function out = computeCellIntegral(obj)
      out = abs(obj.mixedProduct)/2;
    end
    function M   = computeMassMatrix(obj)

      w = ones(size(obj.cells,2)) + eye(size(obj.cells,2));
      d = obj.cellIntegral()/sum(w(:));
      
      i = [obj.cells(:,1); obj.cells(:,1); obj.cells(:,1); ...
           obj.cells(:,2); obj.cells(:,2); obj.cells(:,2); ...
           obj.cells(:,3); obj.cells(:,3); obj.cells(:,3)];
      j = [obj.cells(:,1); obj.cells(:,2); obj.cells(:,3); ...
           obj.cells(:,1); obj.cells(:,2); obj.cells(:,3); ...
           obj.cells(:,1); obj.cells(:,2); obj.cells(:,3)];    
      s = [w(1,1)*d  ; w(1,2)*d  ; w(1,3)*d; ...
           w(2,1)*d  ; w(2,2)*d  ; w(2,3)*d; ...
           w(3,1)*d  ; w(3,2)*d  ; w(3,3)*d];
      
      M = sparse(i, j, s, size(obj.points,1), size(obj.points,1), length(i));
      
    end  
    function out = computePhi(obj)
      
      out = cell(3,2);

      out{1,1} = (obj.points(obj.cells(:,2),2) ...
                 -obj.points(obj.cells(:,3),2))./obj.mixedProduct();
      out{2,1} = (obj.points(obj.cells(:,3),2) ...
                 -obj.points(obj.cells(:,1),2))./obj.mixedProduct();
      out{3,1} = (obj.points(obj.cells(:,1),2) ...
                 -obj.points(obj.cells(:,2),2))./obj.mixedProduct();
      out{1,2} = (obj.points(obj.cells(:,3),1) ...
                 -obj.points(obj.cells(:,2),1))./obj.mixedProduct();
      out{2,2} = (obj.points(obj.cells(:,1),1) ...
                 -obj.points(obj.cells(:,3),1))./obj.mixedProduct();
      out{3,2} = (obj.points(obj.cells(:,2),1) ...
                 -obj.points(obj.cells(:,1),1))./obj.mixedProduct();
      
    end  
    
  end
  
end