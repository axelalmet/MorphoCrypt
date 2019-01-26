function [solMeshNew, SolNew, gammaNew, EbNew, PNew] = UpdateStandardLinearSolidSolution(solMeshOld, solOld, W, parameters, options)

% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);

y0 = parameters.y0;
beta = parameters.beta;
nu = parameters.nu;
K = parameters.K;
dt = parameters.dt;
Eb = parameters.Eb;
g = parameters.g;

% Spring stresses
POld = parameters.P;

% Define new gamma
gammaOld = parameters.gamma;
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

EbNew = Eb;

% Define the ODEs
Odes = @(x, M) StandardLinearSolidFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

solMeshNew = Sol.x;
SolNew = Sol;

% solMeshNew = solMeshOld;

% SolNew.x = solMeshOld;
% SolNew.y = deval(Sol, solMeshOld);

SNew = deval(Sol, solMeshOld, 1);
XNew =  deval(Sol, solMeshOld, 2);
YNew =  deval(Sol, solMeshOld, 3);

DeltaNew = sqrt((XNew - SNew).^2 + (YNew).^2);

if (nu == 0)
    
    PNew = (1 - beta)*K.*(DeltaNew - y0);
    
else
    
    PNew = (1 - dt*beta/nu)^(-1).*(POld + dt*beta*(1 - beta)*K/nu.*(DeltaNew - y0) ...
        + (K/nu).*(DeltaNew - DeltaOld));
    
end
