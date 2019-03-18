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

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

currentArcLength = parameters.currentArcLength;

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

% Update the new boundary
sloughedAmount = trapz(0:dt:t, 0.25*g*(1 + tanh(H*((0:dt:t) - T))));

splitIndex = find(solOld.x == 1, 1);

% Define new functionx
sloughLfunct = @(l) sloughedAmount - (currentArcLength(end) - ...
    interp1(SOld([1:splitIndex, (splitIndex + 2):end]), currentArcLength([1:splitIndex, (splitIndex + 2):end]), l));

LNew = fsolve(sloughLfunct, 0, optimset('Display','off'));





