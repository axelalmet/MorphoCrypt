function [solNew, gammaNew, EbNew, KNew, PxNew, PyNew, ANew] = UpdateAgeingAndPreSloughingSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

W = parameters.W;
sigma = parameters.sigma;
currentArcLength = parameters.currentArcLength;
nu = parameters.nu;
dt = parameters.dt;
% g = parameters.g;
% b1 = parameters.b1;

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;

% Define new gamma
gammaOld = parameters.gamma;
% gammaNew = gammaOld + g*dt;
gammaNew = gammaOld.*(1 + dt.*W(currentArcLength, sigma));
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

% % Update the ages
AOld = parameters.A;
ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
% ANew = 2*AOld + dt - W(currentArcLength, sigma).*AOld;
% 
% b1 = 0.5*(1 - min(ANew)/max(ANew));
% 
% EbNew = 1 - b1.*W(SOld, sigma);
parameters.Eb = EbNew;

% Define the ODEs
Odes = @(x, M) RemodellingFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters);

% Solve using bvp4c.
solNew = bvp4c(Odes, Bcs, solOld, options);
