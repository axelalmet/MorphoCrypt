function EvolveRemodellingFoundationShapeWithContactAlongRegion
%% Load the solutions
outputDirectory = '../../../Solutions/RemodellingFoundation/';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_homoggrowth_w0_0p04_post_constarea';
load([outputDirectory, 'sols_', outputValues, '.mat'], 'contactSols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationContactSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
load([outputDirectory, 'contactMoments_', outputValues,'.mat'], 'momentAtContact') % Foundation stresses
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

%% Let's try contact along a region now
contactSolNew.x = contactSols{end - 2}(1,:);
contactSolNew.y = contactSols{end - 2}(2:end,:);

parameters.Px = foundationContactSols{end - 2}(1,:);
parameters.Py = foundationContactSols{end - 2}(2,:);

% % First update the foundation
XOld = contactSolNew.y(2,:);
YOld = contactSolNew.y(3,:);
nu = parameters.nu;
dt = parameters.dt;

PxOld = parameters.Px;
PyOld = parameters.Py;

parameters.Px = PxOld + nu*dt*(XOld - PxOld);
parameters.Py = PyOld + nu*dt*(YOld - PyOld);

% Split the solution back into two, before contact and after contact.

splitIndex = find(contactSolNew.x == 1, 1);

initInnerSol.x = contactSolNew.x(1:splitIndex);
initInnerSol.y = contactSolNew.y(:, 1:splitIndex);

initOuterSol.x = contactSolNew.x((splitIndex + 1):end) - 1;
initOuterSol.y = contactSolNew.y(1:7, (splitIndex + 1):end);

% Split up the parameter structures as well
innerParameters = parameters;
PxOld = parameters.Px;
innerParameters.Px = PxOld(1:splitIndex);

PyOld = parameters.Py;
innerParameters.Py = PyOld(1:splitIndex);

outerParameters = parameters;

PxOld = parameters.Px;
outerParameters.Px = PxOld((splitIndex + 1):end);

PyOld = parameters.Py;
outerParameters.Py = PyOld((splitIndex + 1):end);
% 
% gammaOld = parameters.gamma;
% g = parameters.g;

% innerParameters.gamma = gammaOld + g*dt;
% outerParameters.gamma = gammaOld + g*dt;
%% Solve the system for contact along a region

% Define the ODEs and BCs
% InnerRegionDerivFun = @(x, M) RemodellingFoundationWithRepulsionInContactRegionOdes(x, M, initInnerSol, innerParameters);
InnerRegionDerivFun = @(x, M) RemodellingFoundationNormalPressureInContactRegionOdes(x, M, initInnerSol, innerParameters);

% Set the boundary conditions
% InnerRegionBcFun = @(Ml, Mr) SelfContactInContactRegionBCs(Ml, Mr, innerParameters);
InnerRegionBcFun = @(Ml, Mr) SelfContactNormalPressureInContactRegionBCs(Ml, Mr, innerParameters);

maxPoints = 1e4;

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');

tic
% Solve the system.
contactRegionInnerSolOld = bvp4c(InnerRegionDerivFun, InnerRegionBcFun, initInnerSol, solOptions);
toc

% Use the inner solution to feed into the outer solution
outerParameters.sc = contactRegionInnerSolOld.y(8, 1);
initOuterSol.y(1, 1) = outerParameters.sc;

% Define the derivative and BC function
% OuterRegionDerivFun = @(x, M) RemodellingFoundationWithRepulsionOutContactRegionOdes(x, M, initOuterSol, outerParameters);
OuterRegionDerivFun = @(x, M) RemodellingFoundationNormalPressureOutContactRegionOdes(x, M, initOuterSol, outerParameters);

OuterRegionBcFun = @(Ml, Mr) SelfContactOutContactRegionBCs(Ml, Mr, outerParameters);

tic
contactRegionOuterSolOld = bvp4c(OuterRegionDerivFun, OuterRegionBcFun, initOuterSol, solOptions);
toc

%% Interpolate the new solutions on the new meshes

% Before-contact solution
innerPxOld = innerParameters.Px;
innerParameters.Px = interp1(initInnerSol.x, innerPxOld, contactRegionInnerSolOld.x);

innerPyOld = innerParameters.Py;
innerParameters.Py = interp1(initInnerSol.x, innerPyOld, contactRegionInnerSolOld.x);

% After-contact solution
outerPxOld = outerParameters.Px;
outerParameters.Px = interp1(initOuterSol.x, outerPxOld, contactRegionOuterSolOld.x);

outerPyOld = outerParameters.Py;
outerParameters.Py = interp1(initOuterSol.x, outerPyOld, contactRegionOuterSolOld.x);
%%

tMax = 5.0;
tMin = newTimes(end);
dt = parameters.dt;
dt = 0.5*dt;
parameters.dt = dt;
regionContactTimes = tMin:dt:tMax;
numSols = length(regionContactTimes);
contactRegionSols = cell(numSols, 1);
foundationContactRegionSols = cell(numSols, 1);

contactRegionSols{1} = [contactRegionInnerSolOld.x, (1 + contactRegionOuterSolOld.x); ...
    contactRegionInnerSolOld.y, [contactRegionOuterSolOld.y; ...
    contactRegionOuterSolOld.y(1,1).*ones(1, length(contactRegionOuterSolOld.x))]];

foundationContactRegionSols{1} = [innerParameters.Px, outerParameters.Px; ...
    innerParameters.Py, outerParameters.Py];

tic
for i = 2:numSols
    
    [contactRegionInnerSolNew, contactRegionOuterSolNew, EbNew, ...
        gammaNew, innerPxNew, innerPyNew, outerPxNew, outerPyNew] = ...
        UpdateRemodellingFoundationShapeWithRegionContact(contactRegionInnerSolOld, innerParameters, ...
        contactRegionOuterSolOld, outerParameters, solOptions);
    
    innerParameters.gamma = gammaNew;
    outerParameters.gamma = gammaNew;
    
    % Update the foundation shapes
    innerParameters.Px = interp1(contactRegionInnerSolOld.x, innerPxNew, contactRegionInnerSolNew.x);
    innerParameters.Py = interp1(contactRegionInnerSolOld.x, innerPyNew, contactRegionInnerSolNew.x);
    
    outerParameters.Px = interp1(contactRegionOuterSolOld.x, outerPxNew, contactRegionOuterSolNew.x);
    outerParameters.Py = interp1(contactRegionOuterSolOld.x, outerPyNew, contactRegionOuterSolNew.x);
    
    % Update the contact point
    outerParameters.sc = contactRegionInnerSolNew.y(8, 1);
    
    contactRegionSols{i} = [contactRegionInnerSolNew.x, (1 + contactRegionOuterSolNew.x); ...
        contactRegionInnerSolNew.y, [contactRegionOuterSolNew.y; ...
        contactRegionOuterSolNew.y(1,1).*ones(1, length(contactRegionOuterSolNew.x))]];
    
    foundationContactRegionSols{i} = [innerParameters.Px, outerParameters.Px; ...
        innerParameters.Py, outerParameters.Py];
    
    % Plot the shape
    figure(1)
    hold on
    plot(contactRegionInnerSolNew.y(2,:), contactRegionInnerSolNew.y(3,:))
    plot(contactRegionOuterSolNew.y(2,:), contactRegionOuterSolNew.y(3,:))
    
    contactRegionInnerSolOld = contactRegionInnerSolNew;
    contactRegionOuterSolOld = contactRegionOuterSolNew;
    
end

toc

%%

for i = 1:5
    
    figure(1)
    hold on
    plot(contactRegionSols{i}(3,:), contactRegionSols{i}(4,:))
    
end

legend('t = 3.3', 't = 3.325', 't = 3.35', 't = 3.375', 't = 3.4')
