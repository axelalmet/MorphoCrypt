function [SolNew, gammaNew, KNew, PxNew, VNew, LNew] = UpdateUnbuckledSloughingSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);

g = parameters.g;
nu = parameters.nu;
dt = parameters.dt;
t = parameters.t;
H = parameters.H;
T = parameters.T;

% Spring stresses
PxOld = parameters.Px;

W = parameters.W;
sigma = parameters.sigma;

% As = parameters.As;

currentArcLength = parameters.currentArcLength;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + g*dt;
% gammaNew = gammaOld.*(1 + dt.*W(currentArcLength, sigma));
parameters.gamma = gammaNew;

% Set new foundation stiffness
KNew = parameters.K;

% Update the foundation shape
PxNew = PxOld + dt*nu.*(XOld - PxOld);

parameters.Px = PxNew;


% Set the new velocity
VNew = cumtrapz(SOld, g.*ones(1, length(SOld)));

% Define the ODEs
Odes = @(x, M, region) UnbuckledRodRemodellingFoundationSloughingOdes(x, M, region, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) UnbuckledRodSloughingBCs(Ml, Mr, parameters);

% Solve using bvp4c.
SolNew = bvp4c(Odes, Bcs, solOld, options);

% Update the ages
% AOld = parameters.A;
% ANew = 2.*AOld + dt - gammaNew./gammaOld.*AOld;

% Update the new boundary

sloughedAmount = trapz(0:dt:t, 0.25.*(1 + tanh(H.*((0:dt:t) - T))));

splitIndex = find(solOld.x == 1, 1);

% Define new functionx
sloughLfunct = @(l) sloughedAmount - (currentArcLength(end) - ...
    interp1(SOld([1:splitIndex, (splitIndex + 2):end]), currentArcLength([1:splitIndex, (splitIndex + 2):end]), l));

LNew = fsolve(sloughLfunct, 0, optimset('Display','off'));






