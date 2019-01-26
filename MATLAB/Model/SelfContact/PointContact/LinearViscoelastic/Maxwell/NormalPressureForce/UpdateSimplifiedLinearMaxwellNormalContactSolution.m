function [solMeshNew, SolNew, gammaNew, EbNew, PNew, uHatNew] = UpdateSimplifiedLinearMaxwellNormalContactSolution(solMeshOld, solOld, W, parameters, options)

% Get the relevant parameters to update growth in time
YOld = solOld.y(3,:);
mOld = solOld.y(7,:);


b1 = parameters.b1;
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
gammaNew = gammaOld.*(1 + dt*(W(solOld.y(1,:), sigma)));
% gammaNew = gammaOld.*(1 + dt*(W(solOld.y(1,:), sigma./gammaOld)));


EbNew = 1 - b1.*W(solOld.y(1,:), sigma./gammaNew);
% EbNew = 1 - b1.*W(solOld.y(1,:), sigma);
parameters.Eb = EbNew;

% Define the ODEs
Odes = @(x, M, region) SimplifiedLinearMaxwellFoundationNormalContactOdes(x, M, region, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) SelfPointContactNormalPressureBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

distanceMatrix = abs(Sol.y(2,2:end) - Sol.y(2,1:(end - 1)));

% Re-mesh the new solution if the inf-norm is too large
if (max(distanceMatrix(:)) > 0.25) 
    
    contactIndex = find(solMeshOld == 1, 1);
    
    solMeshNew = [linspace(0, 1, 2*contactIndex), linspace(1, 2, 2*(length(solMeshOld) - contactIndex))];   
    
    % For sanity
    solMeshNew(1) = 0;
    solMeshNew(end) = 2;
    
else
    
    solMeshNew = solMeshOld;
    
end


SolNew.x = solMeshNew;
SolNew.y = deval(Sol, solMeshNew);

YNew = SolNew.y(3,:);

size(YNew)
if (length(YNew) > length(YOld))
    YNew = InterpolateToNewMesh(solMeshOld, solMeshNew, YNew);
    size(YNew)
end


PNew = YNew - YOld + (1 - dt*nu).*POld; 

uHatNew = uHatOld + etaK*dt.*mOld./Eb;
