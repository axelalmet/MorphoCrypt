function [SolNew, gammaNew, KNew, PNew, VNew, ANew] = UpdateUnbuckledExtensibleRodPreSloughingSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.x;
XOld = solOld.y(1,:);
nOld = solOld.y(2,:);

AOld = parameters.A;

W = parameters.W;
sigma = parameters.sigma;
nu = parameters.nu;
dt = parameters.dt;
mu = parameters.mu;
ns = parameters.ns;

% Update growth
gammaOld = parameters.gamma;
% gammaNew = gammaOld.*(1 + dt*(W(SOld, sigma) + mu.*(nOld - ns).*(nOld < ns)));
gammaNew = gammaOld.*(1 + dt*(W(XOld, sigma) + mu.*(nOld - ns)));
parameters.gamma = gammaNew;

% Set new foundation stiffness
KNew = parameters.K;

% Set the new velocity
VNew = cumtrapz(solOld.x, (gammaNew - gammaOld)./(dt));

% Define the ODEs
Odes = @(x, M) UnbuckledExtensibleRodWithRemodellingFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) UnbuckledExtensibleRodClampedBCs(Ml, Mr, parameters); 
% Bcs = @(Ml, Mr) UnbuckledExtensibleRodSpringBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);
SolNew.x = Sol.x;
SolNew.y = Sol.y;

% Spring stresses
POld = parameters.P;

% Update the foundation shape
PNew = POld + dt*nu.*(XOld - POld);

% Update the age
ANew = 2*AOld + dt - (gammaNew./gammaOld).*AOld;
