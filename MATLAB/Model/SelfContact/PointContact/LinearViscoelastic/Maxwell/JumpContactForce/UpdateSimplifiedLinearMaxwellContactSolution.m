function [SolNew, gammaNew, PNew, uHatNew] = UpdateSimplifiedLinearMaxwellNormalContactSolution(solMesh, solOld, W, parameters, options)

% Get the relevant parameters to update growth in time
YOld = solOld.y(3,:);
mOld = solOld.y(7,:);


sigma = parameters.sigma;
nu = parameters.nu;
etaK = parameters.etaK;
dt = parameters.dt;
Eb = parameters.Eb;

% Spring stresses
POld = parameters.P;
uHatOld = parameters.uHat;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld.*(1 + dt*(W(solOld.x, sigma)));

% Define the ODEs
Odes = @(x, M, region) SimplifiedLinearMaxwellFoundationContactOdes(x, M, region, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) SelfPointContactBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

SolNew.x = solMesh;
SolNew.y = deval(Sol, solMesh);

YNew = InterpolateToNewMesh(solOld.x, SolNew.x, SolNew.y(3,:));

PNew = YNew - YOld + (1 - dt*nu).*POld; 

uHatNew = uHatOld + etaK*dt.*mOld./Eb;