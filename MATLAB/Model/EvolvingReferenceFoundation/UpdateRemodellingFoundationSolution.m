function [SolNew, gammaNew, EbNew, KNew, PxNew, PyNew, uHatNew] = UpdateRemodellingFoundationSolution(solOld, solGuess, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
FOld = solOld.y(4,:);
GOld = solOld.y(5,:);
thetaOld = solOld.y(6,:);
mOld = solOld.y(7,:);

n3Old = FOld.*cos(thetaOld) + GOld.*sin(thetaOld);

% g = parameters.g;
W = parameters.W;
currentArcLength = parameters.currentArcLength;
sigma = parameters.sigma;
nu = parameters.nu;
g = parameters.g;
dt = parameters.dt;
b1 = parameters.b1;
chi = parameters.chi;
Eb = parameters.Eb;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + g*dt;
% gammaNew = gammaOld.*(1 + dt*W(currentArcLength, sigma));
parameters.gamma = gammaNew;

% Set new bending stiffness
% EbNew = 1 - b1.*W(SOld, sigma);
% parameters.Eb = EbNew;
EbNew = parameters.Eb;

% Set new foundation stiffness
KNew = parameters.K;

% Define the ODEs
% Odes = @(x, M) RemodellingFoundationWithHorizontalRepulsionOdes(x, M, solGuess, parameters);
Odes = @(x, M) RemodellingFoundationOdes(x, M, solGuess, parameters);

% Define the BCs
Bcs = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
SolNew = bvp4c(Odes, Bcs, solOld, options);

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;

% Update the foundation shape
PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);

% Update the intrinsic curvature
uHatOld = parameters.uHat;
uHatNew = uHatOld + chi*dt.*mOld./Eb;
