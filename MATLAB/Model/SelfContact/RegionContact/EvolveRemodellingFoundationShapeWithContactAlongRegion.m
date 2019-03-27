function EvolveRemodellingFoundationShapeWithContactAlongRegion
%% Load the solutions
outputDirectory = '../../../Solutions/RemodellingFoundation/';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_homoggrowth_w0_0p04_post_repulsion_Q_1_N_2';
load([outputDirectory, 'sols_', outputValues, '.mat'], 'contactSols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationContactSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
load([outputDirectory, 'contactMoments_', outputValues,'.mat'], 'momentAtContact') % Foundation stresses
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

%% Let's try contact along a region now
contactSolCurrent.x = contactSols{end - 1}(1,:);
contactSolCurrent.y = contactSols{end - 1}(2:end,:);

contactSolPrevious.x = contactSols{end - 2}(1,:);
contactSolPrevious.y = contactSols{end - 2}(2:end,:);

parameters.Px = foundationContactSols{end - 1}(1,:);
parameters.Py = foundationContactSols{end - 1}(2,:);

%%

% % First update the foundation
XOld = contactSolCurrent.y(2,:);
YOld = contactSolCurrent.y(3,:);
nu = parameters.nu;
dt = parameters.dt;

PxOld = parameters.Px;
PyOld = parameters.Py;
 
parameters.Px = PxOld + nu*dt*(XOld - PxOld);
parameters.Py = PyOld + nu*dt*(YOld - PyOld);

% Construct the new solution
initSol = ConstructNewGuessForMultiPointBVPs(contactSolCurrent, contactSolPrevious);
initSol.y(8,:) = initSol.y(8,:);
initSol.y(9,:) = 0.005.*ones(1, length(initSol.x));

%% Solve the system for contact along a region

% Define the ODEs and BCs
regionDerivFun = @(x, M, region) RemodellingFoundationWithRepulsionSelfContactRegionOdes(x, M, region, initSol, parameters);

% Set the boundary conditions
regionBcFun = @(Ml, Mr) SelfContactRegionBCs(Ml, Mr, initSol, parameters);

maxPoints = 1e4;
% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');

tic
contactRegionSolOld = bvp4c(regionDerivFun, regionBcFun, initSol, solOptions);

toc

%%


% Split the solution back into two, before contact and after contact.
splitIndex = find(initSol.x == 1, 1);

initInnerSol.x = initSol.y(1, 1:splitIndex);
initInnerSol.y = initSol.y(1:7, 1:splitIndex);
initInnerSol.parameters = initSol.y(8, 1);

initOuterSol.x = initSol.y(1, (splitIndex + 1):end);
initOuterSol.y = initSol.y(1:7, (splitIndex + 1):end);
initOuterSol.parameters = initSol.y(8, 1);

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
gammaOld = parameters.gamma;
g = parameters.g;

innerParameters.gamma = gammaOld;
outerParameters.gamma = gammaOld;

%%
% Solve the system for contact along a region

% % % Define the ODEs and BCs
% InnerRegionDerivFun = @(x, M, sc) RemodellingFoundationWithRepulsionInContactRegionOdes(x, M, sc, initInnerSol, innerParameters);
% % InnerRegionDerivFun = @(x, M) RemodellingFoundationNormalPressureInContactRegionOdes(x, M, initInnerSol, innerParameters);
% 
% % Set the boundary conditions
% InnerRegionBcFun = @(Ml, Mr, sc) SelfContactInContactRegionBCs(Ml, Mr, sc, innerParameters);
% % InnerRegionBcFun = @(Ml, Mr) SelfContactNormalPressureInContactRegionBCs(Ml, Mr, innerParameters);
% % 
maxPoints = 1e4;
% % 
% % % Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-6,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');
% % 
% tic
% % % Solve the system.
% contactRegionInnerSolOld = bvp4c(InnerRegionDerivFun, InnerRegionBcFun, initInnerSol, solOptions);
% toc

% Use the inner solution to feed into the outer solution

% Define the derivative and BC function
OuterRegionDerivFun = @(x, M, sc) RemodellingFoundationWithRepulsionOutContactRegionOdes(x, M, sc, initOuterSol, outerParameters);
% OuterRegionDerivFun = @(x, M) RemodellingFoundationNormalPressureOutContactRegionOdes(x, M, initOuterSol, outerParameters);

OuterRegionBcFun = @(Ml, Mr, sc) SelfContactOutContactRegionBCs(Ml, Mr, sc, outerParameters);

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
regionContactTimes = tMin:dt:tMax;
numSols = length(regionContactTimes);
contactRegionSols = cell(numSols, 1);
foundationContactRegionSols = cell(numSols, 1);

contactRegionSols{1} = [contactRegionInnerSolOld.x, contactRegionOuterSolOld.x; ...
    contactRegionInnerSolOld.y, contactRegionOuterSolOld.y];

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
    
    % Update the solution structure
    contactRegionSols{i} = [contactRegionInnerSolNew.x, contactRegionOuterSolNew.x; ...
        contactRegionInnerSolNew.y, contactRegionOuterSolNew.y];
    
    foundationContactRegionSols{i} = [innerParameters.Px, outerParameters.Px; ...
        innerParameters.Py, outerParameters.Py];
    
    % Plot the shape
    figure(1)
    hold on
    plot(contactRegionInnerSolNew.y(1,:), contactRegionInnerSolNew.y(2,:))
    plot(contactRegionOuterSolNew.y(1,:), contactRegionOuterSolNew.y(2,:))
    
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
