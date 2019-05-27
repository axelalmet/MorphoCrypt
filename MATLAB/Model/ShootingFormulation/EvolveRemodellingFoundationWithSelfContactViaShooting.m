function EvolveRemodellingFoundationWithSelfContactViaShooting
%% Load the relevant data
% Load the solutions previously computed
outputValues = 'Eb_0p5_sigmaE_0p005_nu_10_K_50_L_1_homoggrowth_w0_0p02';
outputDirectory = '../../Solutions/Shooting/';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
% load([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

% Initialise the current solution, gamma and the stresses

endIndex = 5;
solFromData.x = Sols{end - endIndex}(1,:);
solFromData.y = Sols{end - endIndex}(2:8,:);

% gammaFromData = gammaSols{end - 1}(2,:);
gammaFromData = parameters.gamma - endIndex*(parameters.g)*(parameters.dt);
% PFromData = foundationSols{end - endIndex};

PxFromData = Sols{end - endIndex}(9,:);
PyFromData = Sols{end - endIndex}(10,:);

EbFromData = Sols{end - endIndex}(end, :);

% currentArcLengthFromData = parameters.currentArcLength;

% Set the "Wnt" function
% parameters.W = @(x, sigma) exp(-((x)./sigma).^2);

%% Construct the initial guess for bvp4c

% Initialise the guesses for the contact point
possibleContactPoints = find(abs(solFromData.y(6,:) + 0.5*pi) < 1e-2);
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

% initSol.y = [initSol.y; sCGuess.*ones(1, length(initSol.y(1,:)))];
% initSol.y = [initSol.y; sCGuess.*ones(1, length(initSol.y(1,:))); AGuess; ...
%     pCGuess.*ones(1, length(initSol.y(1,:)))];

initSol.x = [linspace(0, 1, length(solFromData.x(1:contIndex))), ...
    linspace(1, 2, length(solFromData.x(contIndex:end)))];

initPx = PxFromData([1:contIndex,contIndex:end]);

parameters.Px = initPx;

initPy = PyFromData([1:contIndex, contIndex:end]);

parameters.Py = initPy;

initEb = EbFromData([1:contIndex, contIndex:end]);

parameters.Eb = initEb;

firstGamma = gammaFromData;
parameters.gamma = firstGamma;

%% Solve for the first contact solution

solTol = 1e-4;

% Set the tolerances and max. number of mesh points
solOptions = odeset('RelTol', 1e-6,'AbsTol', 1e-6, 'Vectorized', 'On');
previousSol.y = Sols{end - endIndex + 1}(2:(end - 1),:);
% parameters.dt = 0.5*parameters.dt;

initGuess = [initSol.y([3, 4, 7], 1)' ...
            + parameters.dt*(previousSol.y([3, 4, 7], 1)' - initSol.y([3, 4, 7], 1)'), ...
                sCGuess, 100];

maxTries = 100;

tic
[initUnknowns, solResiduals, dtNew, numTries] = ShootForSelfContactAtAPointSolution(initGuess, initSol, parameters, solOptions, solTol, maxTries);

toc

%%


%% Get the new solution and update growth and the foundation accordingly

dtNew = parameters.dt;
solMesh = initSol.x;
parameters.dt = dtNew;
nu = parameters.nu;
g = parameters.g;
gammaOld = parameters.gamma;

% % Update growth
parameters.gamma = gammaOld + g*dtNew;

% Update the foundation
XOld = initSol.y(2,:);
YOld = initSol.y(3,:);

PxOld = parameters.Px;
PyOld = parameters.Py;

parameters.Px = PxOld + dtNew*nu*(XOld - PxOld);
parameters.Py = PyOld + dtNew*nu*(YOld - PyOld);

tic
initContactSol = GetSelfContactAtAPointSolution(initGuess, solMesh, initSol, parameters, solOptions);

toc
