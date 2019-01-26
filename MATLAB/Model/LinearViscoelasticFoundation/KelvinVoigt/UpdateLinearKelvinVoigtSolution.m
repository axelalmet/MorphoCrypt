function [solMeshNew, SolNew, gammaNew, PNew, uHatNew] = UpdateLinearKelvinVoigtSolution(solMeshOld, solOld, W, parameters, options)

% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
mOld = solOld.y(7,:);

g = parameters.g;
nu = parameters.nu;
dt = parameters.dt;
etaK = parameters.etaK;
Eb = parameters.Eb;

uHatOld = parameters.uHat;

% Spring stresses
DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);

% Define new gamma
gammaOld = parameters.gamma;
% gammaNew = gammaOld.*(1 + dt*(W(solOld.x, sigma) + mu.*(n3Old - n3s)));
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Define the ODEs
Odes = @(x, M) LinearKelvinVoigtFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

solMeshNew = solMeshOld;

SolNew.x = solMeshOld;
SolNew.y = deval(Sol, solMeshOld);

SNew = SolNew.y(1,:);
XNew = SolNew.y(2,:);
YNew = SolNew.y(3,:);
DeltaNew = sqrt((XNew - SNew).^2 + (YNew).^2);

PNew = DeltaNew  + (nu/dt).*(DeltaNew - DeltaOld);
uHatNew = uHatOld + etaK*dt.*mOld./Eb;