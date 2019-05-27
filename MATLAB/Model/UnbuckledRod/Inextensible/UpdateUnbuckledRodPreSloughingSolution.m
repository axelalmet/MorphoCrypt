function [SolNew, gammaNew, KNew, PxNew, VNew] = UpdateUnbuckledRodPreSloughingSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
nxOld = solOld.y(4,:);

W = parameters.W;
currentArcLength = parameters.currentArcLength;
sigma = parameters.sigma;
nu = parameters.nu;
g = parameters.g;
dt = parameters.dt;
mu = parameters.mu;

% Spring stresses
PxOld = parameters.Px;

% Update the foundation shape
PxNew = PxOld + dt*nu.*(XOld - PxOld);

parameters.Px = PxNew;

% Update growth
gammaOld = parameters.gamma;
% gammaNew = gammaOld.*(1 + dt*W(currentArcLength, sigma));
% gammaNew = gammaOld*(1 + g*dt);
gammaNew = gammaOld.*(1 + dt*(g + mu*nxOld));
parameters.gamma = gammaNew;

% Set new foundation stiffness
KNew = parameters.K;

% Set the new velocity
% VNew = cumtrapz(SOld, gammaOld.*W(currentArcLength, sigma));
VNew = cumtrapz(SOld, g.*ones(1, length(SOld)));

% Define the ODEs
Odes = @(x, M) UnbuckledRemodellingFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) UnbuckledRodBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);
SolNew.x = 0:1e-2:1;
SolNew.y = deval(Sol, SolNew.x);
