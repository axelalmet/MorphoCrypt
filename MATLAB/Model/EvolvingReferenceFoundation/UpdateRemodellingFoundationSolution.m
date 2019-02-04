function [SolNew, gammaNew, EbNew, KNew, PxNew, PyNew] = UpdateRemodellingFoundationSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
FOld = solOld.y(4,:);
GOld = solOld.y(5,:);
thetaOld = solOld.y(6,:);

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
n3s = parameters.n3s;

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

% gammaNew = gammaOld.*(1 + dt*(W(currentArcLength, sigma) + mu3.*tanh((n3Old - n3s))));
% gammaNew = gammaOld.*(1 + dt*(W(currentArcLength, sigma) + mu3.*(n3Old - n3s)));
% gammaNew = gammaOld.*(1 + dt*(W(currentArcLength, sigma).*(1 + mu3.*tanh(n3Old - n3s))));
% gammaNew = gammaOld.*(1 + dt*(W(currentArcLength, sigma).*(1 + tanh(mu3*(abs(trapz(SOld, n3Old)) - n3s)).*(n3Old - n3s)./max(abs(n3Old) - n3s))));
% gammaNew = gammaOld.*(1 + dt*(W(currentArcLength, sigma).*(1 + tanh(mu3*(abs(trapz(SOld, n3Old)) - n3s)).*(n3Old - n3s)./max(abs(n3Old) - n3s))));
% gammaNew = gammaOld.*(1 + dt*(W(currentArcLength, sigma) + mu3.*(n3Old - n3s)./(1 + (n3Old - n3s).^2)));
gammaNew = gammaOld.*(1 + dt*(W(currentArcLength, sigma)./trapz(SOld, W(currentArcLength, sigma))));
% gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Set new bending stiffness
% EbNew = 1 - b1.*W(SOld, sigma);
% parameters.Eb = EbNew;
EbNew = parameters.Eb;

% Set new foundation stiffness
KNew = parameters.K;

% Define the ODEs
Odes = @(x, M) RemodellingFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
SolNew = bvp4c(Odes, Bcs, solOld, options);
