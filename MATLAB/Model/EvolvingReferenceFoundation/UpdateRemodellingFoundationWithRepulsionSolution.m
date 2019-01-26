function [SolNew, gammaNew, EbNew, KNew, PxNew, PyNew] = UpdateRemodellingFoundationWithRepulsionSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

g = parameters.g;
sigma = parameters.sigma;
nu = parameters.nu;
dt = parameters.dt;
b1 = parameters.b1;
W = parameters.W;

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;

% Update the foundation shape
PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);

parameters.Px = PxNew;
parameters.Py = PyNew;

% Define new gamma
gammaOld = parameters.gamma;

gammaNew = gammaOld + dt*g;
parameters.gamma = gammaNew;

% Set new bending stiffness
% EbNew = 1 - b1.*W(SOld, sigma);
EbNew = parameters.Eb;

% Set new foundation stiffness
KNew = parameters.K;

% Define the ODEs
Odes = @(x, M) RemodellingFoundationWithHorizontalRepulsionOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
SolNew = bvp4c(Odes, Bcs, solOld, options);
