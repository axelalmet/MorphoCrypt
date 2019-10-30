function [solMeshNew, SolNew, PNew] = UpdateStandardLinearSolidSolutionWithoutGrowth(solMeshOld, solOld, parameters, options)

% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);

beta = parameters.beta;
y0 = parameters.y0;
nu = parameters.nu;
K = parameters.K;
dt = parameters.dt;

% Spring stresses
POld = parameters.P;

% Define the ODEs
Odes = @(x, M) StandardLinearSolidFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters);

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

solMeshNew = Sol.x;
SolNew = Sol;

% solMeshNew = solMeshOld;

% SolNew.x = solMeshOld;
% SolNew.y = deval(Sol, solMeshOld);

SNew = deval(Sol, solMeshOld, 1);
XNew =  deval(Sol, solMeshOld, 2);
YNew =  deval(Sol, solMeshOld, 3);

DeltaNew = sqrt((XNew - SNew).^2 + (YNew).^2);

PNew = (1 + dt*beta*nu).*(POld) + dt*beta*(1 - beta)*K*nu.*(DeltaOld - y0) ...
    + (K).*(DeltaNew - DeltaOld);

