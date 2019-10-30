function EvolveRemodellingFoundationShapeWithContact
% Load the solutions previously computed
% outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_homoggrowth_w0_0';
% outputValues = 'Eb_1_nu_10_k_50_L0_0p125_homoggrowth_w0_0p02';
outputValues = 'Eb_1_nu_10_k_50_L0_0p125_homoggrowth_w0_0p02_repulsion_Q_6_N_2';
outputDirectory = '../../../../Solutions/RemodellingFoundation/';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

% Initialise the current solution, gamma and the stresses

endIndex = 4;
solFromData.x = Sols{end - endIndex}(1,:);
solFromData.y = Sols{end - endIndex}(2:8,:);

% gammaFromData = gammaSols{end - 1}(2,:);
gammaFromData = parameters.gamma - endIndex*(parameters.g)*(parameters.dt);
PFromData =  Sols{end - endIndex}(9:10,:);

PxFromData = PFromData(1,:);
PyFromData = PFromData(2,:);

% EbFromData = Sols{end - 1}(end, :);

% currentArcLengthFromData = parameters.currentArcLength;

% Set the "Wnt" function
parameters.W = @(x, sigma) exp(-((x)./sigma).^2);

parameters.w0 = 0.02;

parameters.dt = 0.01;

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

% Define the solution on the new mesh
initSol.y = [solFromData.y(1:7, 1:contIndex), solFromData.y(1:7, contIndex:end)];

initSol.y = [initSol.y; sCGuess.*ones(1, length(initSol.y(1,:)))];
% initSol.y = [initSol.y; sCGuess.*ones(1, length(initSol.y(1,:))); AGuess; ...
%     pCGuess.*ones(1, length(initSol.y(1,:)))];

initSol.x = [linspace(0, 1, length(solFromData.x(1:contIndex))), ...
    linspace(1, 2, length(solFromData.x(contIndex:end)))];

% Re-define the foundation on the new mesh
initPx = PxFromData([1:contIndex,contIndex:end]);
parameters.Px = initPx;

initPy = PyFromData([1:contIndex, contIndex:end]);
parameters.Py = initPy;

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

% Update the foundation
XOld = InterpolateToNewMesh(initSol.x, initSol.y(2,:), contactSolOld.x);
PxOldOld = parameters.Px; % Weird notation but it makes things easier
PxOld = InterpolateToNewMesh(initSol.x, PxOldOld, contactSolOld.x);

nu = parameters.nu;
dt = parameters.dt;
PxNew = PxOld + nu*dt*(XOld - PxOld);
parameters.Px = PxNew;

YOld = InterpolateToNewMesh(initSol.x, initSol.y(3,:), contactSolOld.x);
PyOldOld = parameters.Py;
PyOld = InterpolateToNewMesh(initSol.x, PyOldOld, contactSolOld.x);
PyNew = PyOld + nu*dt*(YOld - PyOld);
parameters.Py = PyNew;

tMax = 20.0;
tMin = times(end - endIndex);
dt = 0.01;
parameters.dt = dt;
newTimes = tMin:dt:tMax;
numSols = length(newTimes);
momentAtContact = zeros(1, numSols);
contactSols = cell(numSols, 1);

contactSols{1} = [initSol.x; initSol.y; PxOldOld; PyOldOld];
contactSols{2} = [contactSolOld.x; contactSolOld.y; parameters.Px; parameters.Py];

momentAtContact(1) = initSol.y(7, find(initSol.x == 1, 1));
momentAtContact(2) = contactSolOld.y(7, find(contactSolOld.x == 1, 1));

tic
for i = 3:numSols
    
    % Construct the solution guess using the previous two known solutions,
    % if we can.
    if (i == 2)
        solGuess = contactSolOld;
    else
        contactSolPrevious.x = contactSols{i - 2}(1,:);
        contactSolPrevious.y = contactSols{i - 2}(2:(end - 2),:);
        
        contactSolCurrent.x = contactSols{i - 1}(1,:);
        contactSolCurrent.y = contactSols{i - 1}(2:(end - 2),:);
        
        solGuess = ConstructNewGuessForMultiPointBVPs(contactSolCurrent, contactSolPrevious);
        
    end
    
    [contactSolNew, EbNew, gammaNew, PxNew, PyNew] = UpdateRemodellingFoundationShapeWithContact(contactSolOld, solGuess, parameters, solOptions);
    
    % Interpolate the bending stiffness to the new mesh
    %     parameters.Eb =  InterpolateToNewMesh(contactSolOld.x, EbNew, contactSolNew.x);
    
    % Update the foundation shape
    parameters.Px = InterpolateToNewMesh(contactSolOld.x, PxNew, contactSolNew.x);
    parameters.Py = InterpolateToNewMesh(contactSolOld.x, PyNew, contactSolNew.x);
    
    contactSols{i} = [contactSolNew.x; contactSolNew.y; parameters.Px; parameters.Py];
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

% Save the solutions
outputDirectory = '../../../../Solutions/RemodellingFoundation/';
% outputValues = 'Eb_1_nu_10_k_50_L0_0p125_homoggrowth_w0_0p02_post_constarea';
outputValues = 'Eb_1_nu_10_k_50_L0_0p125_homoggrowth_w0_0p02_post_repulsion_Q_6_N_2';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'contactSols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
% save([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationContactSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
save([outputDirectory, 'contactMoments_', outputValues,'.mat'], 'momentAtContact') % Foundation stresses
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times