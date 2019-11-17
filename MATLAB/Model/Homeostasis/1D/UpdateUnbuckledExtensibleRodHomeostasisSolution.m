function [SolNew, gammaNew, PNew, VgNew, VcNew] = UpdateUnbuckledExtensibleRodHomeostasisSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
XOld = solOld.y(1,:);

LMu = parameters.LMu;
nu = parameters.nu;
dt = parameters.dt;
t = parameters.t;
g = parameters.g;
T = parameters.T;

% Update growth
gammaOld = parameters.gamma;
gammaNew = gammaOld.*(1 + dt*g(XOld).*(t <= T));
parameters.gamma = gammaNew;

% Set the new growth velocity
VgNew = cumtrapz(solOld.x, (gammaNew - gammaOld)./(dt));

parameters.L = LMu(t);

% Define the ODEs
Odes = @(x, M) UnbuckledExtensibleRodInHomeostasisODEs(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) UnbuckledExtensibleRodHomeostasisClampedBCs(Ml, Mr, parameters);

% Solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);
SolNew.x = solOld.x;
SolNew.y = deval(Sol, solOld.x);

% Set the new current velocity
XNew = SolNew.y(1,:);

VcNew = cumtrapz(solOld.x, (XNew - XOld)./dt);

% Spring stresses
POld = parameters.P;

% Update the foundation shape
PNew = POld + dt*nu.*(XOld - POld);
