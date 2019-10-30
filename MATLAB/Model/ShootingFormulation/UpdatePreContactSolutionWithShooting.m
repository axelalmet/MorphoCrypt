function [solNew, gammaNew, EbNew, KNew, PxNew, PyNew, uHatNew] = UpdatePreContactSolutionWithShooting(solGuess, solOld, solMesh, parameters, options, solTol)
% Obtain the new solution and relevant quantities due to growth  and a
% remodelling foundation.

% Get the relevant parameters to update growth in time
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

W = parameters.W;
nu = parameters.nu;
g = parameters.g;
dt = parameters.dt;
b1 = parameters.b1;
currentArcLength = parameters.currentArcLength;
sigma = parameters.sigma;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Set new bending stiffness
EbNew = 1 - b1.*W(currentArcLength, sigma);
parameters.Eb = EbNew;

% Set new foundation stiffness
KNew = parameters.K;

% Integrate the system

solUnknowns = ShootForPreContactSolution(solGuess, solOld, parameters, options, solTol);

% Evaluate the solution with the right values
solNew = GetPreContactSolution(solUnknowns, solMesh, solOld, parameters, options);

% Update the intrinsic curvature
uHatNew = parameters.uHat;

% Spring stresses
PxOld = parameters.Px;
PyOld = parameters.Py;

% Update the foundation shape
PxNew = PxOld + dt*nu.*(XOld - PxOld);
PyNew = PyOld + dt*nu.*(YOld - PyOld);
