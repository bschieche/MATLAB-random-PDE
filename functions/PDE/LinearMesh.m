classdef LinearMesh < handle
% Stores mesh data and perform mesh-related computations.
%
% This class serves as a super-class for different types of meshes that
% have the following in common:
%
% - All mesh cells are simplices of the same dimension.
% - All mesh points are corner points of the mesh cells.
%
% This is an abstract class. The abstract methods have to implemented in a
% sub-class.
%
% This is a handle class. With a handle class one can implement a mechanism
% to compute certain properties on demand and store the respective data for
% later queries. (Look at the property get methods for more details.)
%
% Sebastian Ullmann
% Numerical Analysis and Scientific Computing
% Technische Universität Darmstadt
% 29/12/2010
% 
% See also LinearMesh/LinearMesh and the sub-classes kardosmesh2d and
% kardosmesh3d.
  
  properties (SetAccess = private)
  % Mesh data.
  
    points        % mesh points
    cells         % mesh cell vertices (indices referring to points)  
  
  end
  
  properties (SetAccess = private, Dependent = true)   
  % Size of the mesh data.
    
    nDim          % number of dimensions
    nPoints       % number of mesh points
    nCells        % number of mesh cells
    nCorners      % number of corners per mesh cell
    
  end
  
  properties (SetAccess = private, Transient = true)
  % Properties that are computed on demand and not stored in files.
    
    cellIntegral     % integral of the unit function over each mesh cell
    massMatrix       % finite element mass matrix
    lumpedMassMatrix % lumped finite element mass matrix
    mixedProduct     % signed cell volume/area/length
    ddx              % spatial gradient
    delta            % maximum cell length
    phi              % spatial derivatives of FEM basis functions
    
  end

  methods
      
    function obj = LinearMesh(points,cells,nDim)
    %LinearMesh Create an instance of LinearMesh. Use via subclass only.
    %
    %Syntax:
    %
    %  obj = LinearMesh(points,cells,nDim)
    %
    %Arguments:
    %
    %  points  An array of mesh points.
    %
    %  cells   An array of triangles. (Indices into points.)
    %
    %  nDim    Number of spatial dimensions.
    %
    %  obj     The LinearMesh object.
    %
    %Remark:
    %
    %  The constructor can only be called from a constructor of a subclass
    %  that specializes dimension-dependent behavior by supplying
    %  implementations for abstract methods.
    %
    %See also LinearMesh and the subclasses kardosmesh2d and kardosmesh3d.
      
      % Check if the number of spatial dimensions is correct.
      assert(size(points,2) == nDim  , ...
        ['size(points,2) should be ' num2str(nDim) ...
         ', but is ' num2str(size(points,2)) '.']);
     
      % Check if the mesh consists of simplices.
      assert(size(cells,2) == nDim+1, ...
        ['size(cells,2) should be ' num2str(nDim+1) ...
         ', but is ' num2str(size(cells,2)) '.']);
     
      % Check if the indices are one-based.
      assert(min(cells(:)) >= 1, ...
        ['min(cells(:)) should be >= 1, but is ' num2str(min(cells(:))) '.']);
    
      % Check if there number of points equals the maximum index.
      assert(max(cells(:)) <= size(points,1), ...
        ['max(cells(:)) should be <= ' num2str(size(points,1)) ...
         ', but is ' num2str(max(cells(:))) '.']);
      
      % Assing the properties.
      obj.points = points;
      obj.cells  = cells;
      
    end
    function out = get.nDim(obj)
      out = size(obj.points,2);
    end
    function out = get.nPoints(obj)
      out = size(obj.points,1);
    end
    function out = get.nCells(obj)
      out = size(obj.cells,1);
    end
    function out = get.nCorners(obj)
      out = size(obj.cells,2);
    end
    function out = get.cellIntegral(obj)
      if isempty(obj.cellIntegral)
        obj.cellIntegral = obj.computeCellIntegral();
      end
      out = obj.cellIntegral;
    end    
    function out = get.massMatrix(obj)
      if isempty(obj.massMatrix)
        obj.massMatrix = obj.computeMassMatrix();
      end
      out = obj.massMatrix;      
    end   
    function out = get.lumpedMassMatrix(obj)
      if isempty(obj.lumpedMassMatrix)
        obj.lumpedMassMatrix = obj.computeLumpedMassMatrix();
      end
      out = obj.lumpedMassMatrix;      
    end   
    function out = get.mixedProduct(obj)
      if isempty(obj.mixedProduct)
        obj.mixedProduct = obj.computeMixedProduct();
      end
      out = obj.mixedProduct;
    end
    function out = get.ddx(obj)
      if isempty(obj.ddx)
        obj.ddx = obj.computeDdx();
      end
      out = obj.ddx;      
    end
    function out = get.delta(obj)
      if isempty(obj.delta)
        obj.delta = obj.computeDelta();
      end
      out = obj.delta;
    end
    function out = get.phi(obj)
      if isempty(obj.phi)
        obj.phi = obj.computePhi();
      end
      out = obj.phi;
    end
    
  end
  
  methods (Sealed = true)
  % The implementations of these methods do not depend on the spatial
  % dimension, so they are defined here.
    
    function out = grad(obj,u)
    %grad Compute gradient of scalar or vector.
      if isnumeric(u)
        assert(size(u,1) == obj.nPoints);
        out = cell(obj.nDim,1);
        for i=1:obj.nDim
          out{i} = obj.ddx{i}*u;
        end
      elseif iscell(u) && numel(u) == 1;
        out = obj.grad(u{:});
      elseif iscell(u) && all(size(u) == [obj.nDim 1])
        out = cell(obj.nDim,obj.nDim);
        for j=1:obj.nDim
          out(:,j) = obj.grad(u{j});
        end
      else
        error('Gradient can not be computed. Check dimensions of input.')
      end
    end
    function out = gradT(obj,u)
    %gradT Compute transpose of gradient of scalar or vector.
      g = obj.grad(u);
      out = cell([size(g,2) size(g,1)]);
      for i=1:size(g,1)
        for j=1:size(g,2)
          out(j,i) = g(i,j);
        end
      end
    end
    function out = div(obj,u)
    %div Compute divergence of vector or matrix.
      assert(iscell(u));
      if all(size(u) == [obj.nDim 1])
        out = zeros(size(obj.cells,1),size(u{1},2));
        for i=1:obj.nDim
          out = out + obj.ddx{i}*u{i};
        end
      elseif all(size(u) == [obj.nDim obj.nDim])
        out = cell(obj.nDim,1);
        for i=1:obj.nDim
          out{i} = obj.div(u(:,i));
        end
      else
        error('Divergence can not be computed. Check dimensions of input.')
      end
    end 
    function [d idx] = findPoint(obj,x)
    %FINDPOINT Find mesh node closest to a given point in space.
      s = zeros(size(obj.points,1),1);
      for i=1:size(obj.points,2)
        s = s+(x(i)-obj.points(:,i)).^2;
      end
      [d idx] = min(s);
    end
    function out = scale(obj)
    %scale Compute the extremal coordinates.
      out = [];
      for i=1:size(obj.points,2)
        out = [out; min(obj.points(:,i)); max(obj.points(:,i))];              %#ok<AGROW>
      end
    end
    function out = L2NormSquared(obj,u)
    %Compute the square of the spatial L2-norm of u.
    
      if isnumeric(u)
        out = sum(u.*(obj.massMatrix*u));
      elseif iscell(u)
        out = zeros(1,size(u{1},2));
        for d=1:numel(u)
          out = out + obj.L2NormSquared(u{d});
        end
      else
        error(['Input of class ' class(u) ' not allowed']);
      end

    end
    function out = L2ErrorSquared(obj,u,uh)
    %Compute the square of the spatial L2-norm of the error e=u-uh.

      if isnumeric(u) && isnumeric(uh)
        out = obj.L2NormSquared(bsxfun(@minus,u,uh));
      elseif iscell(u) && iscell(uh)
        out = zeros(1,max(size(u{1},2),size(uh{1},2)));  
        for d=1:numel(u)
          out = out + obj.L2ErrorSquared(u{d},uh{d});
        end
      else
        error(['Input of class ' class(u) ...
               ' (u) and class ' class(uh) ' (uh) not allowed']);
      end

    end
    
  end
  
  methods (Abstract = true, Access = protected)
  % The implementations of these methods depend on the spatial dimension,
  % so they have to be implemented in sub-classes.
    
    out = computeMassMatrix(obj)
    out = computeCellIntegral(obj)    
    out = computeMixedProduct(obj)
    out = computeDdx(obj)
    out = computePhi(obj)
    
  end
  
  methods (Access = private)
    
    function out = computeDelta(obj)
      out = inf(size(obj.cells,1),1);
      for i=1:size(obj.cells,2)
        pi = obj.points(obj.cells(:,i),:);
        for j=i+1:size(obj.cells,2)
          pj = obj.points(obj.cells(:,j),:);
          out = min(out,sum((pi-pj).^2,2));
        end
      end
      out = sqrt(out);
    end
    function out = computeLumpedMassMatrix(obj)
    %computeLumpedMass Compute the finite element lumped mass matrix.
    %
    %Syntax:
    %
    %  out = computeLumpedMass(obj)
    %
    %Arguments:
    %
    %  out    The finite element lumped mass matrix.
    %
    %  obj    The LinearMesh object.
    
      s = repmat(obj.cellIntegral/size(obj.cells,2), size(obj.cells,2),1);
      out = sparse(obj.cells(:), obj.cells(:), s, size(obj.points,1), size(obj.points,1));
    end
  end
  
  methods (Access = private, Static)
    function varargout = adjustNCols(varargin)
      assert(nargin == nargout);
      maxNCols = 0;
      for i=1:length(varargin)
        if size(varargin{i},2) > maxNCols
          maxNCols = size(varargin{i},2);
        end
      end
      for i=1:length(varargin)
        if size(varargin{i},2) == maxNCols
          varargout{i} = varargin{i};
        elseif size(varargin{i},2) == 1
          varargout{i} = repmat(varargin{i},[1 maxNCols]);
        else
          error('Dimensions do not fit.');
        end
      end
    end
  end
  
  
end