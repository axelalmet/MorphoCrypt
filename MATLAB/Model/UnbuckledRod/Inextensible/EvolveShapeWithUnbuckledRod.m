function EvolveShapeWithUnbuckledRod
% Set the parameters
kf = 0.01; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
L = 1;
K = 12*kf*L0^4/(w*h^3);
dt = 0.05; % Time step

% Get the initial solution from AUTO
solData = load('../../../../Data/planarmorphorodsinextkf0p01L1sol_1');
% solData = load('../../../Data/seashellsK50L1sol_1');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,[2:5])';

% % Load previous heterogeneous solution to obtain a better guess for the solution
% outputDirectory = '../../Solutions/RemodellingFoundation/';
% outputValues = 'Eb_0p5_init_sigmaE_3w_nu_10_kf_0p01_L0_0p125_homoggrowth';
% load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% solFromData.x = Sols{2}(1,:);
% solFromData.y = Sols{2}(2:(end - 1), :);

sigma = 2*w/L0; % "Width" of Wnt gradient
w0 = 0.0*sigma;
% sigma = 0.005;

% % Define the Wnt function
W = @(S, width) exp(-((S)./width).^2);

parameters.W = W;

eta = 1.0/24; % Growth timescale (24 hours)
etaF = eta^(-1);

g = 1;
mu = 0.01;

parameters.g = g;
parameters.currentArcLength = solFromData.y(1,:);
parameters.K = K;% Foundation stiffness
parameters.L = 0.5*L; % Rod length
parameters.sigma = sigma; % Width of wnt gradient
parameters.nu = 0.5*etaF*eta; % Foundation relaxation timescale
parameters.dt = dt; % Time step
parameters.mu = mu;

% Shift the solutions so that we only focus on the 2nd half.
solFromData.x = solFromData.x(300:end) - 0.5;
solFromData.x(1) = 0;
solFromData.x = solFromData.x/(solFromData.x(end));
solFromData.y = solFromData.y(:, 300:end);
solFromData.y(4,:) = solData(1:302,5)';
solFromData.y(1, :) = solFromData.y(1, :) - 0.5;
solFromData.y(2, :) = solFromData.y(2, :) - 0.5;
%


%% Solve the initial bvp to obtain a structure for the first solution.
SOld = solFromData.y(1,:);

% Initialise growth
% firstGamma = 1 + dt*W(solFromData.y(1,:), sigma);
% firstGamma = (1 + dt*g).*(SOld < 0.66*SOld(end));
firstGamma = 1 + g*dt;
parameters.gamma = firstGamma;

parameters.l = trapz(SOld, firstGamma.*ones(1, length(SOld)));

solFromData.y(2,:) = firstGamma.*solFromData.y(1,:);
solFromData.y(3,:) = 0.*solFromData.y(1,:);
solFromData.y(4,:) = 0.*solFromData.y(1,:);

% Initialise foundation shape
parameters.Px = SOld;

% Define the ODEs and BCs
DerivFun = @(x, M) UnbuckledRemodellingFoundationOdes(x, M, solFromData, parameters);

% Set the boundary conditions
BcFun = @(Ml, Mr) UnbuckledRodBCs(Ml, Mr, parameters);

maxPoints = 1e6;

tic
% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');

% Solve the system.
numSol = bvp4c(DerivFun, BcFun, solFromData, solOptions);

toc

initSol = numSol;

%%
% Update the initial data
solOld = initSol;

initS = solOld.y(1,:);

% Interpolate all heterogeneous quantities to the new initial mesh
% parameters.Eb = 1 - b1.*W(initS, sigma);
% parameters.gamma = interp1(solFromData.x, firstGamma, initSol.x);
parameters.currentArcLength = cumtrapz(initS, parameters.gamma.*ones(1, length(initS)));
% parameters.currentArcLength = initS.*(parameters.gamma);

parameters.Px = initS;

parameters.V = cumtrapz(initS, g.*ones(1, length(initS)));

TMax = 1000.0;
times = 0:dt:TMax;
numSols = length(times);

% Initialise the solutions
Sols = cell(numSols, 1);

%  The first solution is always flat
flatSol.x = initSol.x;
flatSol.y = initSol.y;

flatSol.y(3,:) = 0.*flatSol.y(3,:);
flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);

Sols{1} = [flatSol.x; flatSol.y; zeros(1, length(flatSol.x)); zeros(1, length(flatSol.x))];

% First non-trivial solution
Sols{2} = [initSol.x; initSol.y; parameters.gamma.*ones(1, length(initSol.x)); parameters.Px; parameters.currentArcLength; parameters.V];

% Update the solutions in time

tic 
for i = 3:numSols
    
    parameters.t = times(i);
    
%     solCurrent.x = Sols{i - 1}(1,:);
%     solCurrent.y = Sols{i - 1}(2:end,:);
%     
%     solPrevious.x = Sols{i - 2}(1,:);
%     solPrevious.y = Sols{i - 2}(2:end,:);
%     
%     solGuess = ConstructNewGuess(solCurrent, solPrevious);
    
    % Update the solution
    [solNew, gammaNew, KNew, PxNew, VNew] = UpdateUnbuckledRodPreSloughingSolution(solOld, parameters, solOptions);
    
    %     Update the incremental growth
    gammaOld = parameters.gamma;
    gammaInc = (gammaNew - gammaOld)./(dt*gammaOld);
%     gammaIncremental = gammaInc;
    gammaIncremental = interp1(solOld.y(1,:), gammaInc, solNew.y(1,:));
    
    % Update gamma to the new mesh
    parameters.gamma = interp1(solOld.y(1,:), gammaNew, solNew.y(1,:));
%         parameters.gamma = gammaNew;
    
    % Update the foundation shape
    parameters.Px = interp1(solOld.y(1,:), PxNew, solNew.y(1,:));
    
    % Update the current arclength
    parameters.currentArcLength = cumtrapz(solNew.y(1,:), parameters.gamma.*ones(1, length(solNew.x)));
    
    % Update the velocity
    parameters.V = interp1(solOld.x, VNew, solNew.x);
    
    % Stop the solution if the curve self-intersects
    if (times(i) > 1.0 )
        
        Sols = Sols(1:(i - 1));
        times = times(1:(i - 1));
%         foundationSols = foundationSols(1:(i - 1));
        
        break
    end
    
    solOld = solNew;
    
%         Sols{i} = [solOld.x; solOld.y; gammaIncremental];
%     Sols{i} = [solOld.x; solOld.y; parameters.Eb];
    Sols{i} = [solOld.x; solOld.y; gammaIncremental.*ones(1, length(solOld.x)); parameters.Px; parameters.currentArcLength; parameters.V];
    
end

toc

%%
% Sols = Sols(1:(i - 1));
% times = times(1:(i - 1));
% foundationSols = foundationSols(1:(i - 1));

% Save the solutions
outputDirectory = '../../Solutions/UnbuckledRod/';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_homoggrowth_g_0p4_presloughing_T_2p5';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

%% Set up solutions for sloughing

dt = parameters.dt;
g = parameters.g;
LOld = parameters.L;
gammaOld = parameters.gamma;
% AOld = parameters.A;

% gammaNew = gammaOld.*(1 + dt*W(parameters.currentArcLength, parameters.sigma));
gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

% ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
% parameters.A = ANew;

%%

parameters.T = 2.5;
H = 10;

parameters.currentArcLength = cumtrapz(solNew.y(1,:), parameters.gamma.*ones(1, length(solNew.x)));

sloughedAmount = trapz(0:dt:parameters.t, 0.25.*(1 + tanh(H*((0:dt:parameters.t) - parameters.T))));

% Define new functionx
sloughLfunct = @(l) sloughedAmount - (parameters.currentArcLength(end) - interp1(solOld.y(1,:), parameters.currentArcLength, l));

LNew = fsolve(sloughLfunct, 0, optimset('Display','off'));

%%

parameters.L = LNew;

% Initialise the times
TMax = 5.0;
newTimes = (times(end) + dt):dt:TMax;
numSols = length(newTimes);

% Initialise new solution structures
sloughingSols = cell(numSols, 1);
cellSloughingPositions = cell(numSols, 1);

sloughSolOld.x = [solOld.x, 1 2];
sloughSolOld.y = [solOld.y, repmat(solOld.y(:, end), [1 2])];

% Need to modify the solution components individually
sloughSolOld.y(1, 1:(end - 2)) = solOld.y(1, :)./(solOld.y(1, end)).*LNew;
sloughSolOld.y(1, end - 1) = LNew; % Impose new sloughing region

PxOld = parameters.Px;
parameters.Px = [PxOld, repmat(PxOld(end), [1 2])];

parameters.currentArcLength = sloughSolOld.y(1,:).*parameters.gamma;

%% Simulate now with sloughing
tic

trackedPoints = 0:0.1:0.5;

parameters.T = 2.5;
parameters.H = 10;

% Update the solutions in time
for i = 1:numSols
    
    parameters.t = newTimes(i);
    
    % Update the solution
    [sloughSolNew, gammaNew, KNew, PxNew, VNew, LNew] = UpdateUnbuckledSloughingSolution(sloughSolOld, parameters, solOptions);
    
    
    % Update the solutions, gamma, and the foundation shape
%     parameters.gamma = InterpolateToNewMesh(sloughSolOld.x, gammaNew, sloughSolNew.x);
        parameters.gamma = gammaNew;
    parameters.Px = InterpolateToNewMesh(sloughSolOld.x, PxNew, sloughSolNew.x);
    parameters.L = LNew;
%     parameters.A = InterpolateToNewMesh(sloughSolOld.x, ANew, sloughSolNew.x);
%     parameters.A = ANew;

    parameters.currentArcLength = cumtrapz(sloughSolNew.y(1,:), parameters.gamma.*ones(1, length(sloughSolNew.x)));
    
    parameters.V = InterpolateToNewMesh(sloughSolOld.x, VNew, sloughSolNew.x);
    
    sloughSolOld = sloughSolNew;
    
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if ( (~isempty(InterX([sloughSolNew.y(2,:); sloughSolNew.y(3,:)])))) % Stopping criterion
        
        sloughingSols = sloughingSols(1:(i - 1));
        %         gammaSloughingSols = gammaSloughingSols(1:(i - 1));
        newTimes = newTimes(1:(i - 1));
        
        break
    end
    
    sloughingSols{i} = [sloughSolOld.x; sloughSolOld.y; parameters.Px; parameters.V];
    
    splitIndex = find(sloughSolOld.x == 1, 1);
    % Find where the tagged cells are
    taggedPositions = zeros(length(trackedPoints),2);
    for j = 1:length(trackedPoints)
        taggedPositions(j,:) = [interp1(sloughSolOld.y(1, [1:splitIndex, (splitIndex + 2):end]), ...
            sloughSolOld.y(2, [1:splitIndex, (splitIndex + 2):end]), trackedPoints(j)), ...
            interp1(sloughSolOld.y(1, [1:splitIndex, (splitIndex + 2):end]), ...
            sloughSolOld.y(3, [1:splitIndex, (splitIndex + 2):end]), trackedPoints(j))];
    end
    
    cellSloughingPositions{i} = taggedPositions;
    
end

toc

%% Save teh solutions
fullSols = [Sols; sloughingSols];
% fullGammaSols = [Sols; gammaSloughingSols];
fullTimes = [times, newTimes];

fullCellPositions = cellSloughingPositions;

%% Save the solutions
outputDirectory = '../../Solutions/UnbuckledRod/';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_homoggrowth_g_0p4_timestepsloughing_T_2p5';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'fullSols') % Solutions
% save([outputDirectory, 'gamma_', outputValues,'.mat'], 'fullGammaSols') % Gamma
% save([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'fullFoundationSols') % Foundation stresses
save([outputDirectory, 'taggedcells_', outputValues,'.mat'], 'fullCellPositions') % Migrating cells
save([outputDirectory, 'times_', outputValues, '.mat'], 'fullTimes') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

%%  Examine the migration velocities

migrationPositions = zeros(length(trackedPoints), length(newTimes));
migrationVelocities = zeros(1, length(trackedPoints));

for i = 1:numSols
    taggedPositions = cellSloughingPositions{i};
    
    for j = 1:size(taggedPositions, 1)
        
        migrationPositions(j, i) = taggedPositions(j, 1);
        
    end
    
end

for i = 1:length(trackedPoints)
    
    % Get the positions which aren't 'NaN' (happens due to issues with
    % interpolation
    migrationIndices = find(~isnan(migrationPositions(i, :)));
    % Fit the velocities to a linear fit
    [coeffs, resids] = polyfit(newTimes(migrationIndices), migrationPositions(i, migrationIndices), 1);
    migrationVelocities(i) = coeffs(1);
end

save([outputDirectory, 'cellvelocities_', outputValues, '.mat'], 'migrationVelocities') % Times

% Load the Kaur and Potten migration data
KaurPottenVelocities = load('~/Documents/MorphoCrypt/Data/KaurPottenMigrationVelocities.csv');

figure
hold on
plot(trackedPoints*125/10, migrationVelocities*(125*24/(24*10)))
plot(KaurPottenVelocities(:, 1), KaurPottenVelocities(:, 2))
    
legend('Homogeneous growth', 'Kaur and Potten (1986)')

