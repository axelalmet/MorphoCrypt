function [SolNew, EbNew, gammaNew, PNew] = UpdateSLSFoundationShapeWithContact(solOld, solGuess, parameters, options)
% Get the relevant parameters to update growth in time
SOld = solOld.y(1,:);
XOld = solOld.y(2,:);
YOld = solOld.y(3,:);

DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);

beta = parameters.beta;
y0 = parameters.y0;
nu = parameters.nu;
K = parameters.K;
dt = parameters.dt;

% Don't need to udpate these at the moment
EbNew = parameters.Eb;
gammaNew = parameters.gamma;

% Spring stresses
POld = parameters.P;

% Define the ODEs
Odes = @(x, M, region) StandardLinearSolidFoundationContactHalfIntervalOdes(x, M, region, solGuess, parameters);

% Define the BCs
Bcs = @(Ml, Mr) SelfPointContactHalfIntervalBCs(Ml, Mr, parameters); 

% Define solve the ODE system using bvp4c.
Sol = bvp4c(Odes, Bcs, solOld, options);

SolNew = Sol;

SNew = deval(Sol, solOld.x, 1);
XNew =  deval(Sol, solOld.x, 2);
YNew =  deval(Sol, solOld.x, 3);

DeltaNew = sqrt((XNew - SNew).^2 + (YNew).^2);

if (nu == 0)
    
    PNew = (1 - beta)*K.*(DeltaNew - y0);
    
else
    
    PNew = (1 - dt*beta/nu)^(-1).*(POld + dt*beta*(1 - beta)*K/nu.*(DeltaNew - y0) ...
        + (K/nu).*(DeltaNew - DeltaOld));
    
end
