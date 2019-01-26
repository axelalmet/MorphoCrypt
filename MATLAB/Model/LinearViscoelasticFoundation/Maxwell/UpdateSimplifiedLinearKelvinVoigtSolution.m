function [solMeshNew, SolNew, gammaNew, PNew, uHatNew] = UpdateSimplifiedLinearKelvinVoigtSolution(solMeshOld, solOld, W, parameters, options)

% Get the relevant parameters to update growth in time
YOld = solOld.y(3,:);
mOld = solOld.y(7,:);

g = parameters.g;
nu = parameters.nu;
dt = parameters.dt;
etaK = parameters.etaK;
Eb = parameters.Eb;

uHatOld = parameters.uHat;

% Spring stresses
y0 = parameters.y0;

% Define new gamma
gammaOld = parameters.gamma;
% gammaNew = gammaOld.*(1 + dt*(W(solOld.x, sigma) + mu.*(n3Old - n3s)));
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% Define the ODEs
Odes = @(x, M) SimplifiedLinearKelvinVoigtFoundationOdes(x, M, solOld, parameters);

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

YNew = SolNew.y(3,:);

% YOld = interp1(solMeshOld, YOld, solMeshNew);

PNew = YNew - y0  + nu*(dt)^(-1).*(YNew - YOld);
uHatNew = uHatOld + etaK*dt.*mOld./Eb;