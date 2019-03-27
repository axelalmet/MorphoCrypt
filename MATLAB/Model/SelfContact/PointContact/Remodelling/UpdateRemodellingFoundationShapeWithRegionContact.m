function [innerSolNew, outerSolNew, EbNew, ...
            gammaNew, innerPxNew, innerPyNew, ...
            outerPxNew, outerPyNew] = UpdateRemodellingFoundationShapeWithRegionContact(innerSolOld, innerParameters, outerSolOld, outerParameters, options)
% Get the relevant parameters to update growth and the foundation shape in time
g = innerParameters.g;
nu = innerParameters.nu;
dt = innerParameters.dt;

% Update the before-contact foundation
innerXOld = innerSolOld.y(1,:);
innerYOld = innerSolOld.y(2,:);

innerPxOld = innerParameters.Px;
innerPyOld = innerParameters.Py;

innerPxNew = innerPxOld + dt*nu.*(innerXOld - innerPxOld);
innerPyNew = innerPyOld + dt*nu.*(innerYOld - innerPyOld);

innerParameters.Px = innerPxNew;
innerParameters.Py = innerPyNew;

% Update the after-contact foundation
outerXOld = outerSolOld.y(1,:);
outerYOld = outerSolOld.y(2,:);

outerPxOld = outerParameters.Px;
outerPyOld = outerParameters.Py;

outerPxNew = outerPxOld + dt*nu.*(outerXOld - outerPxOld);
outerPyNew = outerPyOld + dt*nu.*(outerYOld - outerPyOld);

outerParameters.Px = outerPxNew;
outerParameters.Py = outerPyNew;

% Define new growth
gammaOld = innerParameters.gamma;
gammaNew = gammaOld + g*dt;
innerParameters.gamma = gammaNew;
outerParameters.gamma = gammaNew;

% Set new bending stiffness
EbNew = innerParameters.Eb;

% First solve the before-contact system
innerSolOld.y(2,:) = innerSolOld.y(2,:) + 0.1;
%Define the ODes
InnerOdes = @(x, M) RemodellingFoundationWithRepulsionInContactRegionOdes(x, M, innerSolOld, innerParameters);

% Set the boundary conditions 
InnerBcs = @(Ml, Mr) SelfContactInContactRegionBCs(Ml, Mr, innerParameters);

tic
% Define solve the ODE system using bvp4c.
innerSolNew = bvp4c(InnerOdes, InnerBcs, innerSolOld, options);
toc

%Define the ODes
OuterOdes = @(x, M) RemodellingFoundationWithRepulsionOutContactRegionOdes(x, M, outerSolOld, outerParameters);

% Set the boundary conditions 
OuterBcs = @(Ml, Mr) SelfContactOutContactRegionBCs(Ml, Mr, outerParameters);

tic
% Define solve the ODE system using bvp4c.
outerSolNew = bvp4c(OuterOdes, OuterBcs, outerSolOld, options);
toc