function [SolNew, EbNew, gammaNew, PxNew, PyNew] = UpdateRemodellingFoundationShapeWithRegionContact(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

g = parameters.g;
currentArcLength = parameters.currentArcLength;
sigma = parameters.sigma;
nu = parameters.nu;
dt = parameters.dt;
b1 = parameters.b1;

% Define the anonymous function for Wnt
W = parameters.W;

PxOld = parameters.Px;
PyOld = parameters.Py;

% Update the foundation
PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);

parameters.Px = PxNew;
parameters.Py = PyNew;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Set new bending stiffness
% EbNew = 1 - b1.*W(SOld, sigma);
% parameters.Eb = EbNew;

EbNew = parameters.Eb;

%Define the ODes
Odes = @(x, M, region) RemodellingFoundationWithRepulsionContactRegionOdes(x, M, region, solOld, parameters);

% Set the boundary conditions 
Bcs = @(Ml, Mr) SelfPointContactRegionBCs(Ml, Mr, parameters);

% Define solve the ODE system using bvp4c.
SolNew = bvp4c(Odes, Bcs, solOld, options);