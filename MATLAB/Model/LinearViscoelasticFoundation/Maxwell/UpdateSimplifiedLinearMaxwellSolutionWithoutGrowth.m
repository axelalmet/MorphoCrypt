function [solMeshNew, SolNew, PNew, uHatNew] = UpdateSimplifiedLinearMaxwellSolutionWithoutGrowth(solMeshOld, solOld, parameters, options)

% Get the relevant parameters to update growth in time
YOld = solOld.y(3,:);
mOld = solOld.y(7,:);

nu = parameters.nu;
etaK = parameters.etaK;
dt = parameters.dt;
Eb = parameters.Eb;

% Spring stresses
POld = parameters.P;
uHatOld = parameters.uHat;

% Define the ODEs
Odes = @(x, M) SimplifiedLinearMaxwellFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

distanceMatrix = abs(Sol.y(2:3,2:end) - Sol.y(2:3,1:(end-1)));

% Re-mesh the new solution if the inf-norm is too large
if (max(distanceMatrix(:)) > norm(abs(solMeshOld(2:end) - solMeshOld(1:(end - 1))),inf))    
    solMeshNew = linspace(0, 1, 2*length(solMeshOld));   
    
    % For sanity
    solMeshNew(1) = 0;
    solMeshNew(end) = 1;
else
    solMeshNew = solMeshOld;
end

SolNew.x = solMeshNew;
SolNew.y = deval(Sol, solMeshNew);

YNew = SolNew.y(3,:);

PNew = YNew - YOld + (1 - dt*nu).*POld; 
uHatNew = uHatOld + etaK*dt.*mOld./Eb;