function [SolNew, PxNew, PyNew, uHatNew] = UpdateRemodellingFoundationSolutionWithoutGrowth(solOld, parameters, options)
% Get the relevant parameters to update growth in time
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
mOld = solOld.y(7,:);

Eb = parameters.Eb;
chi = parameters.chi;
nu = parameters.nu;
dt = parameters.dt;

% Define the ODEs
Odes = @(x, M) RemodellingFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
SolNew = bvp4c(Odes, Bcs, solOld, options);

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;

PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);

% Intrinsic curvature
uHatOld = parameters.uHat;
uHatNew = uHatOld + chi*dt.*mOld./Eb;
