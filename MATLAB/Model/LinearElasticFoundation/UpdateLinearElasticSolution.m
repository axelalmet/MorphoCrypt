function [solMeshNew, SolNew, gammaNew, uHatNew] = UpdateLinearElasticSolution(solMeshOld, solOld, W, parameters, options)

% Get the previous solution
mOld = solOld.y(7,:);

% Get the relevant parameters
sigma = parameters.sigma;
Eb = parameters.Eb;
etaK = parameters.etaK;
dt = parameters.dt;

uHatOld = parameters.uHat;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + dt;
% gammaNew = gammaOld.*(1 + dt*(W(solOld.y(1,:), sigma));

parameters.gamma = gammaNew;

% Define the ODEs
Odes = @(x, M) LinearElasticFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 
% Bcs = @(Ml, Mr) HalfIntervalBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

solMeshNew = Sol.x;
SolNew = Sol;

uHatNew = uHatOld + etaK*dt.*mOld./Eb;
