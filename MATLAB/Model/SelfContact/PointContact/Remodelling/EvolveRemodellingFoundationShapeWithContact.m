function EvolveRemodellingFoundationShapeWithContact
% Load the solutions previously computed
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_homoggrowth_w0_0p04_repulsion_Q_10_N_2';
outputDirectory = '../../../../Solutions/RemodellingFoundation/';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

% Initialise the current solution, gamma and the stresses

endIndex = 5;
solFromData.x = Sols{end - endIndex}(1,:);
solFromData.y = Sols{end - endIndex}(2:(end - 1),:);

% gammaFromData = gammaSols{end - 1}(2,:);
gammaFromData = parameters.gamma - endIndex*(parameters.g)*(parameters.dt);
PFromData = foundationSols{end - endIndex};

PxFromData = PFromData(1,:);
PyFromData = PFromData(2,:);

% EbFromData = Sols{end - 1}(end, :);

% currentArcLengthFromData = parameters.currentArcLength;

% Set the "Wnt" function
parameters.W = @(x, sigma) exp(-((x)./sigma).^2);

%% Construct the initial guess for bvp4c

% Initialise the guesses for the contact point
possibleContactPoints = find(abs(solFromData.y(6,:) + 0.5*pi) < 5*1e-2);
contIndex = possibleContactPoints(end);
sCGuess = solFromData.y(1, contIndex);
pCGuess = 0;
AGuess = [cumtrapz(solFromData.y(1, 1:contIndex), ...
    (parameters.gamma).*(-solFromData.y(2, 1:contIndex)).*sin(solFromData.y(6, 1:contIndex))), ...
    trapz(solFromData.y(1, 1:contIndex), ...
    (parameters.gamma).*(-solFromData.y(2, 1:contIndex)).*sin(solFromData.y(6, 1:contIndex)))...
    .*ones(1, length(contIndex:length(solFromData.x)))];

A0 = AGuess(contIndex);

parameters.A0 = A0;

initSol.y = [solFromData.y(1:7, 1:contIndex), solFromData.y(1:7, contIndex:end)];

initSol.y = [initSol.y; sCGuess.*ones(1, length(initSol.y(1,:)))];
% initSol.y = [initSol.y; sCGuess.*ones(1, length(initSol.y(1,:))); AGuess; ...
%     pCGuess.*ones(1, length(initSol.y(1,:)))];

initSol.x = [linspace(0, 1, length(solFromData.x(1:contIndex))), ...
    linspace(1, 2, length(solFromData.x(contIndex:end)))];

initX = initSol.y(2,:);
initPx = PxFromData([1:contIndex,contIndex:end]);

parameters.Px = initPx + (parameters.dt)*(parameters.nu).*(initX - initPx);

initY = initSol.y(3,:);
initPy = PyFromData([1:contIndex, contIndex:end]);

parameters.Py = initPy + (parameters.dt)*(parameters.nu).*(initY - initPy);

% initEb = EbFromData([1:contIndex, contIndex:end]);

% parameters.Eb = initEb;
%
firstGamma = gammaFromData;
parameters.gamma = firstGamma + (parameters.g)*(parameters.dt);

% parameters.currentArcLength = parameters.gamma.*initSol.y(1,:);

%%

% Solve the new contact system

% Define the ODEs and BCs
% DerivFun = @(x, M, region) RemodellingFoundationNormalPressureContactHalfIntervalOdes(x, M, region, initSol, parameters);
DerivFun = @(x, M, region) RemodellingFoundationContactWithRepulsionOdes(x, M, region, initSol, parameters);
% DerivFun = @(x, M, region) RemodellingFoundationContactHalfIntervalOdes(x, M, region, initSol, parameters);

% Set the boundary conditions
% BcFun = @(Ml, Mr) SelfPointContactNormalPressureBCs(Ml, Mr, parameters);
BcFun = @(Ml, Mr) SelfPointContactHalfIntervalBCs(Ml, Mr, parameters);

maxPoints = 1e4;

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');

tic
% Solve the system.
contactSolOld = bvp4c(DerivFun, BcFun, initSol, solOptions);

toc

%% Interpolate the new solutions on the new meshes

% EbOld = parameters.Eb;
% parameters.Eb = InterpolateToNewMesh(initSol.x, EbOld, contactSolOld.x);

PxOld = parameters.Px;
parameters.Px = InterpolateToNewMesh(initSol.x, PxOld, contactSolOld.x);

PyOld = parameters.Py;
parameters.Py = InterpolateToNewMesh(initSol.x, PyOld, contactSolOld.x);

tMax = 5.0;
tMin = times(end - endIndex);
dt = 0.05;
parameters.dt = dt;
newTimes = tMin:dt:tMax;
numSols = length(newTimes);
momentAtContact = zeros(1, numSols);
contactSols = cell(numSols, 1);
foundationContactSols = cell(numSols, 1);

contactSols{1} = [contactSolOld.x; contactSolOld.y];

foundationContactSols{1} = [parameters.Px; parameters.Py];

momentAtContact(1) = contactSolOld.y(7, find(contactSolOld.x == 1, 1));

tic
for i = 2:numSols
    
    % Construct the solution guess using the previous two known solutions,
    % if we can.
    if (i == 2)
        solGuess = contactSolOld;
    else
        contactSolPrevious.x = contactSols{i - 2}(1,:);
        contactSolPrevious.y = contactSols{i - 2}(2:end,:);
        
        contactSolCurrent.x = contactSols{i - 1}(1,:);
        contactSolCurrent.y = contactSols{i - 1}(2:end,:);
        
        solGuess = ConstructNewGuessForMultiPointBVPs(contactSolCurrent, contactSolPrevious);
        
    end
    
    [contactSolNew, EbNew, gammaNew, PxNew, PyNew] = UpdateRemodellingFoundationShapeWithContact(contactSolOld, solGuess, parameters, solOptions);
    
    % Interpolate the bending stiffness to the new mesh
    %     parameters.Eb =  InterpolateToNewMesh(contactSolOld.x, EbNew, contactSolNew.x);
    
    % Update the foundation shape
    parameters.Px = InterpolateToNewMesh(contactSolOld.x, PxNew, contactSolNew.x);
    parameters.Py = InterpolateToNewMesh(contactSolOld.x, PyNew, contactSolNew.x);
    
    contactSols{i} = [contactSolNew.x; contactSolNew.y];
    foundationContactSols{i} = [parameters.Px; parameters.Py];
    momentAtContact(i) = contactSolNew.y(7, find(contactSolNew.x == 1, 1));
    
    %     % Plot the shape
    %     figure(1)
    %     subplot(1, 2, 1)
    %     hold on
    %     plot(contactSolNew.y(2,:), -contactSolNew.y(3,:))
    %
    %     % Plot the bending moment (curvature)
    %     figure(1)
    %     subplot(1, 2, 2)
    %     hold on
    %     plot(contactSolNew.y(1,:), contactSolNew.y(7,:))
    
    if (HasRodHitRegionContact(contactSolNew, parameters) )
        
        contactSols = contactSols(1:(i));
        newTimes = newTimes(1:(i));
        momentAtContact = momentAtContact(1:(i));
        foundationContactSols = foundationContactSols(1:(i));
        
        break
        
    end
    
    contactSolOld = contactSolNew;
    parameters.gamma = gammaNew;
    
end

toc

% %%
% contactSols = contactSols(1:(end - 1));
% newTimes = newTimes(1:(end - 1));
% momentAtContact = momentAtContact(1:(end - 1));
% foundationContactSols = foundationContactSols(1:(end - 1));

%% Save the solutions
outputDirectory = '../../../../Solutions/RemodellingFoundation/';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_homoggrowth_w0_0p04_post_repulsion_Q_10_N_2';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'contactSols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
save([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationContactSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
save([outputDirectory, 'contactMoments_', outputValues,'.mat'], 'momentAtContact') % Foundation stresses
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

%% Load the solutions
outputDirectory = '../../../../Solutions/RemodellingFoundation/';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_homoggrowth_w0_0p04_post_repulsion_Q_1_N_2';
load([outputDirectory, 'sols_', outputValues, '.mat'], 'contactSols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationContactSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
load([outputDirectory, 'contactMoments_', outputValues,'.mat'], 'momentAtContact') % Foundation stresses
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
%%
figure(3)
hold on
% plot(solFromData.y(2,:), -solFromData.y(3,:), 'b','linewidth', 1.5)
% plot(-solFromData.y(2,:), -solFromData.y(3,:), 'b','linewidth', 1.5)
plot(contactSols{1}(3,:), -contactSols{2}(4,:), 'r', 'linewidth', 1.5)
plot(-contactSols{1}(3,:), -contactSols{2}(4,:), 'r', 'linewidth', 1.5)
plot(contactSols{end}(3,:), -contactSols{end}(4,:), 'k', 'linewidth', 1.5)
plot(-contactSols{end}(3,:), -contactSols{end}(4,:), 'k', 'linewidth', 1.5)

%% Let's try contact along a region now
contactSolNew.x = contactSols{end}(1,:);
contactSolNew.y = contactSols{end}(2:end,:);

% First update the foundation
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
initOuterSol.y = contactSolNew.y(:, (splitIndex + 1):end);

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

gammaOld = parameters.gamma;
g = parameters.g;

innerParameters.gamma = gammaOld + g*dt;
outerParameters.gamma = gammaOld + g*dt;
%% Solve the system for contact along a region

% Define the ODEs and BCs
InnerRegionDerivFun = @(x, M) RemodellingFoundationWithRepulsionInContactRegionOdes(x, M, initInnerSol, innerParameters);

% Set the boundary conditions
InnerRegionBcFun = @(Ml, Mr) SelfContactInContactRegionBCs(Ml, Mr, innerParameters);

maxPoints = 1e4;

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');

tic
% Solve the system.
contactRegionInnerSolOld = bvp4c(InnerRegionDerivFun, InnerRegionBcFun, initInnerSol, solOptions);
toc

% Use the inner solution to feed into the outer solution
outerParameters.sc = contactRegionInnerSolOld.y(8, 1);
initOuterSol.y = initOuterSol.y(1:7,:);
initOuterSol.y(1, 1) = outerParameters.sc;

% Define the derivative and BC function
OuterRegionDerivFun = @(x, M) RemodellingFoundationWithRepulsionOutContactRegionOdes(x, M, initOuterSol, outerParameters);
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
