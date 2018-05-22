close all

addpath(genpath('functions'));

problemfolder = '.'; % name of the solution folder
refinementlevel = 3;
sigma = 0.2;

% generate the folder if it doesn't exist
if ~exist([problemfolder,'/solutions'],'file')
    mkdir([problemfolder,'/'],'solutions');
end

sol = readSol(2,1,'finiteElementMesh.geo'); % read coarse grid
points = sol.points(2:3,:)'; % extract the points
cells = sol.triangles(2:4,:)'+1; % extract the triangles

L = 4;
for i = 1:refinementlevel
    [points,cells,~,I] = RefineTriangularGrid(points,cells);
end
mesh = LinearTriangularMesh(points,cells); % construct a mesh with various information
eigInfoName = ['eigInfoL',int2str(L),'refine',int2str(refinementlevel),'.mat'];
if ~exist(eigInfoName,'file')
    % run the eigenvalue compuation by means of a Galerkin approach (on the same meth
    [V,lambda] = GalerkinEigenProblem2(L,points,cells,'full');
    save(eigInfoName,'V','lambda');
else
    load(eigInfoName,'V','lambda'); % load KL data
end

% dimension of random space (number of eigenvalues which explain at least 95% of variance
if length(lambda) < size(points,1) % in case lambda had been truncated before
    d = length(lambda);
else
    d = find(cumsum(lambda)/sum(lambda)>0.95,1,'first');
    lambda = lambda(1:d);
    V = V(:,1:d);
end

for i = 1:d
    figure
    trimesh(cells,points(:,1),points(:,2),V(:,i),'FaceColor','interp','EdgeColor','interp');
    view(2)
    axis equal
    axis tight
end

V = bsxfun(@times,V,sigma*sqrt(lambda')); % multiply the eigenvectors with the square root of the eigenvalues
transformationFunc = @(x)[ones(size(x,1),1) x]; % simple transformation needed for the simulation
tspan = [0 1]; % time horizon
type = 2; % type = 1 for steady state problem
odeSolver = @ode15s; % MATLAB ODE solver for stiff problems
odeOptions = odeset('RelTol',1e-5, 'AbsTol',1e-5);
% generate object for the problem under consideration
heatObj = HeatClass(mesh,[ones(mesh.nPoints,1) V],transformationFunc,tspan,odeSolver,odeOptions,type);

tol = 1.0e-3;
% cell array of function handles for the quantities of interest
qoiHandles = {@(x)qoiSpaceTimeIntegral(x,heatObj.mesh,heatObj.freeNodes); ...
    @(x)qoiTimeIntegral(x,heatObj.mesh,heatObj.freeNodes,0,0); ...
    @(x)qoiSpaceIntegral(x,heatObj.mesh,heatObj.freeNodes); ...
    @(x)qoiMaxValue(x,heatObj.mesh,heatObj.freeNodes,0,0);
    @(x)qoiMaxValue(x,heatObj.mesh,heatObj.freeNodes,0,0.5);
    @(x)qoiMaxValue(x,heatObj.mesh,heatObj.freeNodes,1,0.5);
    @(x)qoiMaxValue(x,heatObj.mesh,heatObj.freeNodes,2,0.5);
    @(x)qoiMaxValue(x,heatObj.mesh,heatObj.freeNodes,2,0);
    @(x)qoiMaxValue(x,heatObj.mesh,heatObj.freeNodes,1,-0.6)};

a = -sqrt(3);
b = sqrt(3);
slevel = 6;
rule = {'GaussPattersonLegendre',ones(1,d)};
func = @(nodes,problemfolder)...
    parabolic2D(nodes,problemfolder,@(x)heatObj.compute(x),qoiHandles);
omega = 0.5;
qoi = 2;
[allnodes,level,surplus,intval,nodesol,alpha,deltagrid] = ...
    anisoHierLagrange(a,b,slevel,rule,d,func,tol,omega,problemfolder,100,qoi);
weights = intval'*alpha;
fprintf('\nAprroximation = %1.6e\n',sqrt(intval'*surplus{2}(:,2)-(intval'*surplus{1}(:,2)).^2));
fprintf('\n Solution with respect to 6 spatial reference points:\n\n')
fprintf('  P = %2d  ',length(allnodes));
fprintf('\n');
for i = 4:9
    fprintf(' %1.6f ',(intval'*surplus{1}(:,i))');
    fprintf('\n');
end

lw = 1;
fs = 22;

mypath2 = [problemfolder,'/solutions'];
load([mypath2,'/nodelist']);

d = size(allnodes,2);

t = tspan(1):0.02:tspan(2);
allmean = zeros(length(t),6);
allmean2 = allmean;
lastmean = 0;
lastmean2 = 0;

u = cell(6,1);

for i = 1:size(allnodes,1)
    [c,j] = tolismember(allnodes(i,:),nodelist(:,1:d),10^-8);
    if c
        I(i) = nodelist(j,end);
    end
    sol = load([mypath2,'/solution',int2str(I(i))]);
    
    u = zeros(mesh.nPoints,length(sol.x));
    u(heatObj.freeNodes,:) = sol.y;
    y1 = mesh.valueAt(0,0,u);
    y2 = mesh.valueAt(0,0.5,u);
    y3 = mesh.valueAt(1,0.5,u);
    y4 = mesh.valueAt(2,0.5,u);
    y5 = mesh.valueAt(2,0,u);
    y6 = mesh.valueAt(1,-0.6,u);
    lastmean = lastmean + weights(i)*u(:,end);
    lastmean2 = lastmean2 + weights(i)*u(:,end).^2;
    
    u = zeros(mesh.nPoints,size(sol.idata.dif3d,2)*size(sol.idata.dif3d,3));
    u(heatObj.freeNodes,:) = reshape(sol.idata.dif3d,size(sol.idata.dif3d,1),size(sol.idata.dif3d,2)*size(sol.idata.dif3d,3));
    
    sol.y = y1;
    sol.idata.dif3d = reshape(mesh.valueAt(0,0,u),1,size(sol.idata.dif3d,2),size(sol.idata.dif3d,3));
    allmean(:,1) = allmean(:,1)+ weights(i)*deval(t,sol)';
    allmean2(:,1) = allmean2(:,1)+ weights(i)*(deval(t,sol)').^2;
    
    sol.y = y2;
    sol.idata.dif3d = reshape(mesh.valueAt(0,0.5,u),1,size(sol.idata.dif3d,2),size(sol.idata.dif3d,3));
    allmean(:,2) = allmean(:,2)+ weights(i)*deval(t,sol)';
    allmean2(:,2) = allmean2(:,2)+ weights(i)*(deval(t,sol)').^2;
    
    sol.y = y3;
    sol.idata.dif3d = reshape(mesh.valueAt(1,0.5,u),1,size(sol.idata.dif3d,2),size(sol.idata.dif3d,3));
    allmean(:,3) = allmean(:,3)+ weights(i)*deval(t,sol)';
    allmean2(:,3) = allmean2(:,3)+ weights(i)*(deval(t,sol)').^2;
    
    sol.y = y4;
    sol.idata.dif3d = reshape(mesh.valueAt(2,0.5,u),1,size(sol.idata.dif3d,2),size(sol.idata.dif3d,3));
    allmean(:,4) = allmean(:,4)+ weights(i)*deval(t,sol)';
    allmean2(:,4) = allmean2(:,4)+ weights(i)*(deval(t,sol)').^2;
    
    sol.y = y5;
    sol.idata.dif3d = reshape(mesh.valueAt(2,0,u),1,size(sol.idata.dif3d,2),size(sol.idata.dif3d,3));
    allmean(:,5) = allmean(:,5)+ weights(i)*deval(t,sol)';
    allmean2(:,5) = allmean2(:,5)+ weights(i)*(deval(t,sol)').^2;
    
    sol.y = y6;
    sol.idata.dif3d = reshape(mesh.valueAt(1,-0.6,u),1,size(sol.idata.dif3d,2),size(sol.idata.dif3d,3));
    allmean(:,6) = allmean(:,6)+ weights(i)*deval(t,sol)';
    allmean2(:,6) = allmean2(:,6)+ weights(i)*(deval(t,sol)').^2;
    
end

figure
[~,I] = sort(allmean(end,:),'descend');
errorbar(repmat(t',1,6),allmean(:,I),sqrt(allmean2(:,I)-allmean(:,I).^2),'LineWidth',lw)
axis tight
set(gca,'Fontsize',fs,'LineWidth',lw);
xlabel('time t')
ylabel('solution at spatial reference points')
legendstr = {'(0,0)','(0,0.5)','(1,0.5)','(2,0.5)','(2,0)','(1,-0.6)'};
legendstr = legendstr(I);
legend(legendstr,'Location','NorthWest')

figure
subplot(2,1,1)
trimesh(cells,points(:,1),points(:,2),lastmean,'FaceColor','interp','EdgeColor','interp');
title('expected value','FontSize',fs)
CM = get(gcf,'Colormap');
view(2);axis equal tight off
subplot(2,1,2)
trimesh(cells,points(:,1),points(:,2),sqrt(lastmean2-lastmean.^2),'FaceColor','interp','EdgeColor','interp');
title('standard deviation','FontSize',fs)
colormap(CM);
view(2);axis equal tight off

dimchecker = [1 2];
figure
plotinterp([a a],[b b],allnodes,level,surplus{1}(:,qoi),[],'hier',dimchecker);
hold on
plot(allnodes(:,dimchecker(1)),allnodes(:,dimchecker(2)),'k.','MarkerSize',15)
set(gca,'FontSize',15,'LineWidth',1)
xlabel('y_1');
ylabel('y_2');

densX = rand(10^5,d)*(b-a)+a;
dens = plotinterp([],[],allnodes,level,surplus{1}(:,qoi),densX,'hier',[]);
[densdata,z] = ksdensity(dens,'function','pdf');
figure;
plot(z,densdata,'k','LineWidth',1);
set(gca,'FontSize',15,'LineWidth',1)
set(gca,'FontSize',18);
xlabel('range of possible realizations')
