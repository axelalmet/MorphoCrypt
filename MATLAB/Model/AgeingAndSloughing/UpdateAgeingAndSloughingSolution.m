function [solNew, gammaNew, EbNew, KNew, PxNew, PyNew, ANew, LNew] = UpdateAgeingAndSloughingSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

g = parameters.g;
nu = parameters.nu;
dt = parameters.dt;
t = parameters.t;
H = parameters.H;
T = parameters.T;
LOld = parameters.L;
W = parameters.W;
currentArcLength = parameters.currentArcLength;
sigma = parameters.sigma;

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Set new bending stiffness
EbNew = parameters.Eb;

% Set new foundation stiffness
KNew = parameters.K;

% Update the foundation shape
PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);

parameters.Px = PxNew;
parameters.Py = PyNew;

% Define the ODEs
Odes = @(x, M, region) RemodellingFoundationWithSloughingOdes(x, M, region, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) SloughingBCs(Ml, Mr, parameters);

% Solve using bvp4c.
solNew = bvp4c(Odes, Bcs, solOld, options);

% Update the ages
AOld = parameters.A;
ANew = 2.*AOld + dt - gammaNew./gammaOld.*AOld;

LNew = LOld + dt*g/gammaOld*(-0.25*(1 + tanh(H*(t - T))) + (0.5 - LOld));
parameters.L = LNew;


