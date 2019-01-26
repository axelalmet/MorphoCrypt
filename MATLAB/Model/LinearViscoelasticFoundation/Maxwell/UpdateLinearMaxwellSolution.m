function [solMeshNew, SolNew, gammaNew, PNew, uHatNew] = UpdateLinearMaxwellSolution(solMeshOld, solOld, W, parameters, options)

% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
mOld = solOld.y(7,:);

DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);

nu = parameters.nu;
etaK = parameters.etaK;
dt = parameters.dt;
Eb = parameters.Eb;
g = parameters.g;

% Spring stresses
POld = parameters.P;
uHatOld = parameters.uHat;

% Define new gamma
gammaOld = parameters.gamma;
% gammaNew = gammaOld.*(1 + dt*(W(solOld.y(1,:), sigma) + mu.*(n3Old - n3s)));
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Define the ODEs
Odes = @(x, M) LinearMaxwellFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);
% 
% solMeshNew = Sol.x;
% SolNew = Sol;
solMeshNew = solMeshOld;

SolNew.x = solMeshOld;
SolNew.y = deval(Sol, solMeshOld);
% 
SNew = deval(Sol, solMeshOld, 1);
XNew =  deval(Sol, solMeshOld, 2);
YNew =  deval(Sol, solMeshOld, 3);

DeltaNew = sqrt((XNew - SNew).^2 + (YNew).^2);

PNew = DeltaNew - DeltaOld + (1 - dt*nu).*POld; 

uHatNew = uHatOld + etaK*dt.*mOld./Eb;