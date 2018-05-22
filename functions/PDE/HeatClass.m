classdef HeatClass

% based on code by 
% Sebastian Ullmann, 2010
% Numerical Analysis and Scientific Computing
% Technische Universität Darmstadt

    properties
        mesh;
        M;
        S;
        F;
        adjointF;
        usedNodes;
        freeNodes;
        dirichletNodes;
        cauchyNodes;
		cavityLeft;
		cavityRight;
		cavityMiddle;
        transformationFunc;
        tspan;
        odeSolver;
        odeOptions;
    end
    
    methods
        function obj = HeatClass(mesh,eigenFunctions,transformationFunc,tspan,odeSolver,odeOptions,type)
            
            obj.transformationFunc = transformationFunc;
            obj.tspan = tspan;
            obj.odeSolver = odeSolver;
            obj.odeOptions = odeOptions;
            obj.mesh = mesh;
            
            %obj.dirichletNodes = obj.mesh.boundary;
            %obj.dirichletNodes = unique(obj.dirichletNodes(:));
            bound = obj.mesh.boundary;
            bound = unique(bound(:));
            tol = 1.0e-14;
            left = bound(abs(mesh.points(bound,1)+3)<tol);
            right = bound(abs(mesh.points(bound,1)-3)<tol);
            lower = bound(abs(mesh.points(bound,2)+0.6)<tol);
            upper = bound(abs(mesh.points(bound,2)-1)<tol);
            obj.cavityLeft = bound(abs(mesh.points(bound,1)+1)<tol & mesh.points(bound,2)<=0);
            obj.cavityRight = bound(abs(mesh.points(bound,1)-1)<tol & mesh.points(bound,2)<=0);
            obj.cavityMiddle = bound(abs(mesh.points(bound,2))<tol & mesh.points(bound,1)<=1 & mesh.points(bound,1)>=-1);
            cavity = unique([obj.cavityLeft;obj.cavityRight;obj.cavityMiddle]);
            
            if type == 1
                obj.dirichletNodes = upper;
            elseif type == 2
                obj.dirichletNodes = [];
            end
            
            % distinguish between types of nodes
            obj.usedNodes      = unique(obj.mesh.cells(:));                  % at any triangle
            obj.freeNodes      = setdiff(obj.usedNodes,obj.dirichletNodes);  % degrees of freedom
            obj.cauchyNodes    = cavity;          % Cauchy boundary conditions
            
            % create matrix
            obj.M = obj.mesh.massMatrix(obj.freeNodes,obj.freeNodes);
            obj.S = obj.stiffnessMatrix(eigenFunctions);
            
            %obj.F = obj.volumeForce;
            obj.F = obj.boundaryFlux(cavity);
            obj.adjointF = obj.adjointRHS;
            
        end
        
        function y = compute(obj,nodes)
            
            nodes = obj.transformationFunc(nodes);
            initialSolution = zeros(size(obj.freeNodes));
            
            y = cell(size(nodes,1),1);
            
            for i=1:size(nodes,1)
                
                Si = 0*obj.S{1};
                
                for j=1:size(nodes,2)
                    Si = Si + obj.S{j}*nodes(i,j);
                end
                
                newOdeOptions = odeset('Mass',obj.M, ...
                    'Jacobian',-Si, ...
                    'MassSingular','no');
                y{i} = obj.odeSolver(@(t,y)(-Si*y + obj.F), ...
                    obj.tspan, ...
                    initialSolution, ...
                    odeset(obj.odeOptions,newOdeOptions));
                y{i}.extdata = [];
                
            end
            
        end
        
        function y = adjoint(obj,nodes)
            
            nodes = obj.transformationFunc(nodes);
            
            initialSolution = zeros(size(obj.adjointF));
            y = cell(size(nodes,1),1);
            
            for i=1:size(nodes,1)
                
                Si = 0*obj.S{1};
                
                for j=1:size(nodes,2)
                    Si = Si + obj.S{j}*nodes(i,j);
                end
                newOdeOptions = odeset('Mass',obj.M, ...
                    'Jacobian',-Si, ...
                    'MassSingular','no');
                y{i} = obj.odeSolver(@(t,y)(-Si*y + obj.adjointF), ...
                    obj.tspan, ...
                    initialSolution, ...
                    odeset(obj.odeOptions,newOdeOptions));
                y{i}.extdata = [];
                
            end
        end
        
        function [outx,outy,elements Delta] = getElementPoints(obj)
            
            coordinates = obj.mesh.points;
            elements    = obj.mesh.cells;
            
            x1 = coordinates(elements(:,1),1);
            y1 = coordinates(elements(:,1),2);
            x2 = coordinates(elements(:,2),1);
            y2 = coordinates(elements(:,2),2);
            x3 = coordinates(elements(:,3),1);
            y3 = coordinates(elements(:,3),2);
            
            outx = [x1 x2 x3];
            outy = [y1 y2 y3];
            
            Delta = (x2-x1).*(y3-y1)-(y2-y1).*(x3-x1);
            
        end
        
        function out = reconstruct(obj,y)
            
            coordinates = obj.mesh.points;
            elements    = obj.mesh.cells;
            
            yAll = zeros(size(coordinates,1),1);
            yAll(obj.freeNodes) = y;
            yAll = yAll(elements);
            
            x1 = coordinates(elements(:,1),1);
            y1 = coordinates(elements(:,1),2);
            x2 = coordinates(elements(:,2),1);
            y2 = coordinates(elements(:,2),2);
            x3 = coordinates(elements(:,3),1);
            y3 = coordinates(elements(:,3),2);
            
            Delta = (x2-x1).*(y3-y1)-(y2-y1).*(x3-x1);
            
            a1 = (y2-y3)./Delta; % Gradient of 1st basis function in point 1 of the triangle w.r.t. x
            a2 = (y3-y1)./Delta;
            a3 = (y1-y2)./Delta;
            b1 = (x3-x2)./Delta;
            b2 = (x1-x3)./Delta;
            b3 = (x2-x1)./Delta;         
            
            %dx = (yAll(:,1).*a1+yAll(:,2).*a2+yAll(:,3).*a3).*abs(Delta)*0.5;
            %dy = (yAll(:,1).*b1+yAll(:,2).*b2+yAll(:,3).*b3).*abs(Delta)*0.5;
            
            dx = (yAll(:,1).*a1+yAll(:,2).*a2+yAll(:,3).*a3);
            dy = (yAll(:,1).*b1+yAll(:,2).*b2+yAll(:,3).*b3);
            
                
%             dx = accumarray(elements(:),[dx;dx;dx],[size(coordinates,1) 1])./...
%                  accumarray(elements(:),0.5*abs([Delta;Delta;Delta]),[size(coordinates,1) 1]);
%             
%             dy = accumarray(elements(:),[dy;dy;dy],[size(coordinates,1) 1])./...
%                  accumarray(elements(:),0.5*abs([Delta;Delta;Delta]),[size(coordinates,1) 1]);
%              
            dx = accumarray(elements(:),[dx;dx;dx],[size(coordinates,1) 1])./...
                 accumarray(elements(:),ones(size([dx;dx;dx])),[size(coordinates,1) 1]);
            
            dy = accumarray(elements(:),[dy;dy;dy],[size(coordinates,1) 1])./...
                 accumarray(elements(:),ones(size([dy;dy;dy])),[size(coordinates,1) 1]);
                      
            out = evalReconstruct([x1 x2 x3],[y1 y2 y3],yAll,dx(elements),dy(elements));           
            
        end
        
        function S = getNodeS(obj,nodes)
            
            nodes = obj.transformationFunc(nodes);
            
            initialSolution = zeros(size(obj.freeNodes));
            S = cell(size(nodes,1),1);
            
            for i=1:size(nodes,1)
                
                Si = 0*obj.S{1};
                
                for j=1:size(nodes,2)
                    Si = Si + obj.S{j}*nodes(i,j);
                end
                S{i} = Si;
            end
        end
        
        function M = stiffnessMatrix(obj,f)
            
            M = cell(size(f,2),1);
            
            coordinates = obj.mesh.points;
            elements    = obj.mesh.cells;
            
            fmid = (f(elements(:,1),:) ...
                +f(elements(:,2),:) ...
                +f(elements(:,3),:))/3;
            
            ind1 = [elements(:,1); elements(:,1); elements(:,1); ...
                elements(:,2); elements(:,2); elements(:,2); ...
                elements(:,3); elements(:,3); elements(:,3)];
            
            ind2 = [elements(:,1); elements(:,2); elements(:,3); ...
                elements(:,1); elements(:,2); elements(:,3); ...
                elements(:,1); elements(:,2); elements(:,3)];
            
            x1 = coordinates(elements(:,1),1);
            y1 = coordinates(elements(:,1),2);
            x2 = coordinates(elements(:,2),1);
            y2 = coordinates(elements(:,2),2);
            x3 = coordinates(elements(:,3),1);
            y3 = coordinates(elements(:,3),2);
            
            Delta = (x2-x1).*(y3-y1)-(y2-y1).*(x3-x1);
            
            a1 = (y2-y3)./Delta;
            a2 = (y3-y1)./Delta;
            a3 = (y1-y2)./Delta;
            b1 = (x3-x2)./Delta;
            b2 = (x1-x3)./Delta;
            b3 = (x2-x1)./Delta;
            
            for i=1:numel(M)
                
                m11 = 0.5*abs(Delta).*(a1.*a1+b1.*b1).*fmid(:,i);
                m12 = 0.5*abs(Delta).*(a1.*a2+b1.*b2).*fmid(:,i);
                m13 = 0.5*abs(Delta).*(a1.*a3+b1.*b3).*fmid(:,i);
                m22 = 0.5*abs(Delta).*(a2.*a2+b2.*b2).*fmid(:,i);
                m23 = 0.5*abs(Delta).*(a2.*a3+b2.*b3).*fmid(:,i);
                m33 = 0.5*abs(Delta).*(a3.*a3+b3.*b3).*fmid(:,i);
                
                m = [m11; m12; m13; ...
                    m12; m22; m23; ...
                    m13; m23; m33];
                
                M{i} = sparse(ind1,ind2,m,size(coordinates,1),size(coordinates,1));
                M{i} = M{i}(obj.freeNodes,obj.freeNodes);
                
            end
            
        end
        
        function f = volumeForce(obj)
            
            coordinates = obj.mesh.points;
            elements    = obj.mesh.cells;
            
            f = ones(size(elements,1),1);
            
            x1 = coordinates(elements(:,1),1);
            y1 = coordinates(elements(:,1),2);
            x2 = coordinates(elements(:,2),1);
            y2 = coordinates(elements(:,2),2);
            x3 = coordinates(elements(:,3),1);
            y3 = coordinates(elements(:,3),2);
            
            Delta = (x2-x1).*(y3-y1)-(y2-y1).*(x3-x1);
            
            i = [elements(:,1); elements(:,2); elements(:,3)];
            f = [abs(Delta).*f; abs(Delta).*f; abs(Delta).*f]/6;
            
            f = accumarray(i,f,[size(coordinates,1) 1]);
            
            f = f(obj.freeNodes);
            
        end
        
        function f = boundaryFlux(obj,cavity)
            
            coordinates = obj.mesh.points;
            bound = obj.mesh.boundary;
            bound = bound(ismember(bound(:,1),cavity) & ismember(bound(:,2),cavity),:);
            
            x1 = coordinates(bound(:,1),1);
            y1 = coordinates(bound(:,1),2);
            x2 = coordinates(bound(:,2),1);
            y2 = coordinates(bound(:,2),2);
            
            edgeIntegral = sqrt((x2-x1).^2+(y2-y1).^2)*0.5;
            
            f = sparse(bound,ones(size(bound)),repmat(edgeIntegral,1,2),size(coordinates,1),1);
            
            f = f(obj.freeNodes);
        end
        
        function f = adjointRHS(obj)
            
            coordinates = obj.mesh.points;
            cells = obj.mesh.cells;
            I = find(coordinates(:,1)==0 & coordinates(:,2)==0);
            f = zeros(size(coordinates,1),1);
            f(I) = 1;
            f = f(obj.freeNodes);
        end
    end
    
end

