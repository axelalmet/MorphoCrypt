function [SolNew, gammaNew, KNew, PNew, VNew] = UpdateUnbuckledExtensibleRodPreSloughingSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
XOld = solOld.y(1,:);

W = parameters.W;
sigma = parameters.sigma;
nu = parameters.nu;
g = parameters.g;
dt = parameters.dt;


% Update growth
gammaOld = parameters.gamma;
gammaNew = gammaOld.*(1 + dt*W(solOld.y(1,:), sigma));
% gammaNew = gammaOld*(1 + g*dt);
% gammaNew = gammaOld.*(1 + dt*(g + mu*nxOld));
% gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Set new foundation stiffness
KNew = parameters.K;

% Set the new velocity
VNew = cumtrapz(solOld.x, gammaOld.*W(XOld, sigma));

% Define the ODEs
Odes = @(x, M) UnbuckledExtensibleRodWithRemodellingFoundationOdes(x, M, solOld, parameters);

% Define the BCs
% Bcs = @(Ml, Mr) UnbuckledExtensibleRodClampedBCs(Ml, Mr, parameters); 
Bcs = @(Ml, Mr) UnbuckledExtensibleRodSpringBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);
SolNew.x = Sol.x;
SolNew.y = Sol.y;

% Spring stresses
POld = parameters.P;

% Update the foundation shape
PNew = POld + dt*nu.*(XOld - POld);
