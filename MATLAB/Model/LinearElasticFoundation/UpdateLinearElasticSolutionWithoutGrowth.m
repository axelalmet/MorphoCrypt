function [SolNew, uHatNew] = UpdateLinearElasticSolutionWithoutGrowth(solOld, parameters, options)

% Get the previous solution
mOld = solOld.y(7,:);

% Get the relevant parameters
dt = parameters.dt;
Eb = parameters.Eb;
chi = parameters.chi;
uHatOld = parameters.uHat;

% Update the intrinsic curvature
uHatNew = uHatOld + chi*dt.*mOld./Eb;
parameters.uHat = uHatNew;

% Define the ODEs
Odes = @(x, M) LinearElasticFoundationOdes(x, M, solOld, parameters);

% Define the BCs
% Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 
Bcs = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

SolNew = Sol;

