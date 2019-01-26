function [SolNew, PxNew, PyNew] = UpdateRemodellingFoundationSolutionWithoutGrowth(solOld, parameters, options)
% Get the relevant parameters to update growth in time
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

nu = parameters.nu;
dt = parameters.dt;

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;

% Define the ODEs
Odes = @(x, M) RemodellingFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
SolNew = bvp4c(Odes, Bcs, solOld, options);

% Obtain the solution at the new mesh
% SolNew.x = solMeshNew;
% SolNew.y = deval(Sol, solMeshNew);

PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);