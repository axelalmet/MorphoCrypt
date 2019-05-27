function [SolNew, gammaNew, KNew, PxNew, VNew] = UpdateUnbuckledRodSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);

W = parameters.W;
currentArcLength = parameters.currentArcLength;
sigma = parameters.sigma;
nu = parameters.nu;
g = parameters.g;
dt = parameters.dt;

% Spring stresses
PxOld = parameters.Px;

% Update the foundation shape
PxNew = PxOld + dt*nu.*(XOld - PxOld);

parameters.Px = PxNew;

% Update growth
gammaOld = parameters.gamma;
% gammaNew = gammaOld.*(1 + dt*W(currentArcLength, sigma));
gammaNew = gammaOld + g*dt;
% gammaNew = gammaOld*(1 + g*dt);
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
SolNew = bvp4c(Odes, Bcs, solOld, options);
