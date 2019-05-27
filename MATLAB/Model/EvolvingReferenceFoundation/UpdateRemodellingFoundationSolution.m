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
mu3 = parameters.mu3;
chi = parameters.chi;
Eb = parameters.Eb;

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;

% Update the foundation shape
PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);

parameters.Px = PxNew;
parameters.Py = PyNew;

% Update the intrinsic curvature
uHatOld = parameters.uHat;
uHatNew = uHatOld + chi*dt.*mOld./Eb;
parameters.uHat = uHatNew;

% Define new gamma
% t = parameters.t;
gammaOld = parameters.gamma;

% gammaNew = gammaOld.*(1 + dt*((1 - alpha)*W(currentArcLength, sigma) + 2*alpha.*(W(currentArcLength + 1.175*sigma, sigma) + W(currentArcLength - 1.175*sigma, sigma))));
% gammaNew = gammaOld.*(1 + dt*((1 - alpha)*W(currentArcLength, sigma) + alpha*mu3*W(currentArcLength, sigma).*(n3Old - n3s)));
% gammaNew = gammaOld.*(1 + dt*(W(currentArcLength, sigma)./trapz(SOld, W(currentArcLength, sigma)).*trapz(SOld, 2*(W(currentArcLength + 1.175*sigma, sigma) + W(currentArcLength - 1.175*sigma, sigma)))));
% gammaNew = gammaOld + g*dt;
gammaNew = gammaOld*(1 + g*dt);
parameters.gamma = gammaNew;

% Set new bending stiffness
% EbNew = 1 - b1.*W(currentArcLength, sigma);
% EbNew = 1 - b1.*W(SOld, sigma);
% parameters.Eb = EbNew;
EbNew = parameters.Eb;

% Set new foundation stiffness
KNew = parameters.K;

% Define the ODEs
Odes = @(x, M) RemodellingFoundationOdes(x, M, solGuess, parameters);

% Define the BCs
Bcs = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
SolNew = bvp4c(Odes, Bcs, solOld, options);
