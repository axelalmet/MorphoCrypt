function [SolNew, EbNew, gammaNew, PxNew, PyNew] = UpdateRemodellingFoundationShapeWithContact(solOld, solGuess, parameters, options)
% Get the relevant parameters to update growth in time

g = parameters.g;
dt = parameters.dt;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Set new bending stiffness
EbNew = parameters.Eb;

%Define the ODes
Odes = @(x, M, region) RemodellingFoundationNormalPressureContactHalfIntervalOdes(x, M, region, solOld, parameters);
% Odes = @(x, M, region) RemodellingFoundationContactHalfIntervalOdes(x, M, region, solOld, parameters);
% Odes = @(x, M, region) RemodellingFoundationContactWithRepulsionOdes(x, M, region, solOld, parameters);

% Set the boundary conditions 
Bcs = @(Ml, Mr) SelfPointContactNormalPressureBCs(Ml, Mr, parameters);
% Bcs = @(Ml, Mr) SelfPointContactHalfIntervalBCs(Ml, Mr, parameters);

% Define solve the ODE system using bvp4c.
SolNew = bvp4c(Odes, Bcs, solGuess, options);

% Update the foundation attachments
nu = parameters.nu;
dt = parameters.dt;

XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

PxOld = parameters.Px;
PyOld = parameters.Py;

PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);