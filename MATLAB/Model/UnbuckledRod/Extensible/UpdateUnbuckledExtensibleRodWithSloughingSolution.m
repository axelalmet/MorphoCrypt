function [SolNew, gammaNew, LNew, KNew, PNew, VNew, ANew] = UpdateUnbuckledExtensibleRodWithSloughingSolution(solOld, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
nOld = solOld.y(3,:);

% Get the relevant parmaeters
% W = parameters.W;
% sigma = parameters.sigma;
nu = parameters.nu;
dt = parameters.dt;
% mu = parameters.mu;
% ns = parameters.ns;
% eta = parameters.eta;
% As = parameters.As;
g = parameters.g;

% Update growth
gammaOld = parameters.gamma;
% gammaNew = gammaOld.*(1 + dt*(W(XOld, sigma) + mu.*(nOld - ns).*(nOld < ns)));
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Set new foundation stiffness
KNew = parameters.K;

% Set the new velocity
VNew = cumtrapz(solOld.x, (gammaNew - gammaOld)./(dt).*ones(1, length(solOld.x)));

% Define the ODEs
Odes = @(x, M) UnbuckledExtensibleRodWithRemodellingFoundationAndSloughingOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) UnbuckledExtensibleRodClampedBCsWithSloughing(Ml, Mr, parameters);

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);
SolNew.x = Sol.x;
SolNew.y = Sol.y;

% Spring stresses
POld = parameters.P;

% Update the foundation shape
PNew = POld + dt*nu.*(XOld - POld);

% Update the age
AOld = parameters.A;

ANew = 2*AOld + dt - (gammaNew./gammaOld).*AOld;

% Update the new length
LOld = parameters.L;

% If Age crosses the critical age, we slide the boundary towards the minimum
% age at this threshold
% LMin = parameters.LMin;
% if (LMin < 1)
% if (AOld(end) > As)
%     LNew = LOld - dt*eta*(AOld(end) - As);
% %     LNew = LMin;
%
%     % We need to now slide all quantities across to the new boundary
%     S0New = SOld./SOld(end).*LNew; % Re-adjust new S0
%
%     % We now interpolate the solutions onto this new mesh, with the adjusted L
%     ANewOld = ANew; % Poor terminology but wil make sense
%     ANew = interp1(SOld, ANewOld, S0New); %Update the age
%
%     gammaNewOld = gammaNew;
%     gammaNew = interp1(SOld, gammaNewOld, S0New); % Update the growth
%
%     PNewOld = PNew;
%     PNew = interp1(SOld, PNewOld, S0New); % Foundation
%
%     VNewOld = VNew;
%     VNew = interp1(SOld, VNewOld, S0New); % Velocity
%
% else
%     LNew = LOld;
% end

% Update the new length
t = parameters.t;
T = parameters.T;
if (t > T)
    mu = trapz(0:dt:t, g.*ones(1, length(0:dt:t)).*((0:dt:t) > T));
    LNew = 1 - mu./gammaOld;
else
    LNew = parameters.L;
end
