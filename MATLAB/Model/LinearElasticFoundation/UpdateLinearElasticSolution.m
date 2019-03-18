function [SolNew, gammaNew, uHatNew] = UpdateLinearElasticSolution(solOld, parameters, options)

% Get the previous solution
mOld = solOld.y(7,:);

% Get the relevant parameters
% W = parameters.W;
% sigma = parameters.sigma;
Eb = parameters.Eb;
chi = parameters.chi;
dt = parameters.dt;
g = parameters.g;

uHatOld = parameters.uHat;

% Update the intrinsic curvature
uHatNew = uHatOld + chi*dt.*mOld./Eb;
parameters.uHat = uHatNew;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + g*dt;
% gammaNew = gammaOld.*(1 + dt*(W(solOld.y(1,:), sigma));

parameters.gamma = gammaNew;

% Define the ODEs
Odes = @(x, M) LinearElasticFoundationOdes(x, M, solOld, parameters);

% Define the BCs
% Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 
Bcs = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

SolNew = Sol;

