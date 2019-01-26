function [SolNew, PNew, uHatNew] = UpdateLinearMaxwellSolutionWithoutGrowth(solMesh, solOld, parameters, options)

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

% Spring stresses
POld = parameters.P;
uHatOld = parameters.uHat;

% Define the ODEs
Odes = @(x, M) LinearMaxwellFoundationOdes(x, M, solOld, parameters);

% Define the BCs
Bcs = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

SolNew.x = solMesh;
SolNew.y = deval(Sol, solMesh);

SNew = SolNew.y(1,:);
XNew = SolNew.y(2,:);
YNew = SolNew.y(3,:);

DeltaNew = sqrt((XNew - SNew).^2 + (YNew).^2);

PNew = DeltaNew - DeltaOld + (1 - dt*nu).*POld; 
uHatNew = uHatOld + etaK*dt.*mOld./Eb;