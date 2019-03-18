function EvolveShapeWithAgeingAndSloughing
% Set the parameters
kf = 0.01; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
% L = 2*sqrt(3)*L0/h; %Dimensionless length
L = 1;
y0 = 0*L;
% K = kf*h/(12*w); % Dimensionless foundation stiffness
K = 12*kf*L0^4/(w*h^3);
% K = 50;
Es = 1; % Stretching stiffness
b1 = 0.0; % Bending stiffness
dt = 0.05; % Time step

A = 0;
H = 10;

% Get the initial solution from AUTO
solData = load('../../../Data/planarmorphorodsinextkf0p01L1sol_1');
% solData = load('../../../Data/seashellsK50L1sol_1');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

solFromData.y(3,:) = y0 + solFromData.y(3,:);

sigma = 2*w/L0; % "Width" of Wnt gradient
% sigma = 0.005;

% % Define the Wnt function
W = @(S, width) exp(-((S)./width).^2);
parameters.W = W;

eta = 1.0/24; % Growth timescale (24 hours)
etaF = eta^(-1);

mu3 = 0.1;

g = 1;

T = 2.5;

parameters.mu3 = mu3;
parameters.g = g;
parameters.currentArcLength = solFromData.y(1,:);
parameters.K = K;% Foundation stiffness
parameters.L = 0.5*L; % Rod length
parameters.y0 = y0; % Rod centreline
parameters.sigma = sigma; % Width of wnt gradient
parameters.eta = eta; % Rate of chemical change
parameters.Es = Es; % Stretch stiffness
parameters.b1 = b1;
parameters.Eb = 1; % Bending stiffness
parameters.ext = 0; % Exstensibility
parameters.nu = 10*etaF*eta; % Foundation relaxation timescale
parameters.etan3 = 0*eta*etaF; % Relaxation of target stress
parameters.etaK = 0*eta*etaF; % Curvature relaxation timescale
parameters.dt = dt; % Time step

parameters.t = dt;
parameters.H = H;
parameters.A = A;
parameters.As = 5.0;
parameters.T = T;

% Calculate the critical buckling stress, so that the difference in stress
% quantifies the relative change
PMin = Get1stCriticalStressValue(parameters);
lowerSearchBound = PMin + 1; upperSearchBound = PMin + 100;
searchTol = 1e-4;
[lowerBound, upperBound] = GetIntervalToFind2ndCriticalStress(lowerSearchBound:searchTol:(upperSearchBound), parameters, searchTol);

PCritical = Get2ndCriticalStressValue(parameters, [lowerBound, upperBound]);
n3s = -PCritical;
parameters.n3s = n3s;

% Shift the solutions so that we only focus on the 2nd half.
solFromData.x = solFromData.x(300:end) - 0.5;
solFromData.x(1) = 0;
solFromData.x = solFromData.x/(solFromData.x(end));
solFromData.y = solFromData.y(:, 300:end);
solFromData.y(4,:) = solData(1:302,5)';
solFromData.y(1, :) = solFromData.y(1, :) - 0.5;
solFromData.y(2, :) = solFromData.y(2, :) - 0.5;

solFromData.y([1; 2], 1) = [0; 0];
%

solMesh = 0:1e-3:1;
trackedPoints = [0.1, 0.2, 0.4];

%% Solve the initial bvp to obtain a structure for the first solution.
SOld = solFromData.y(1,:);

% Initialise growth
% firstGamma = 1 + dt*g;
firstGamma = 1 + dt*W(SOld, sigma);
parameters.gamma = firstGamma;

parameters.Eb = 1;

% Initialise foundation shape
parameters.Px = SOld;
parameters.Py = y0.*ones(1, length(solFromData.x));

% Define the ODEs and BCs
DerivFun = @(x, M) RemodellingFoundationOdes(x, M, solFromData, parameters);

% Set the boundary conditions
BcFun = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters);

maxPoints = 1e6;

tic
% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');

% Solve the system.
numSol = bvp4c(DerivFun, BcFun, solFromData, solOptions);

toc

initSol = numSol;

%% Simulate up to the point of sloughing

% Update the initial data
solOld = initSol;

% % Find where the tagged cells are
% taggedPositions = zeros(length(trackedPoints),2);
% for i = 1:length(trackedPoints)
%     solAtMesh = deval(initSol, solMesh, [1:3]);
%     taggedPositions(i,:) = solAtMesh(2:end, find(solAtMesh(1,:) == trackedPoints(i), 1))';
% end

initS = solOld.y(1,:);
initY = solOld.y(3,:);

% gammaOld = interp1(solFromData.x, parameters.gamma, initSol.x);
AOld = dt;
parameters.A = AOld;
parameters.t = 2*dt;

gammaOld = parameters.gamma;
parameters.currentArcLength = cumtrapz(initS, gammaOld);

parameters.Px = initS;
parameters.Py = dt*(parameters.nu).*initY;

% Set the times we want to solve the problem for
TMax = 5.0;
times = 0:dt:TMax;
numSols = length(times);

% Initialise the solutions
Sols = cell(numSols, 1);
gammaSols = cell(numSols, 1);
foundationSols = cell(numSols, 1);
cellPositions = cell(numSols, 1);

%  The first solution is always flat
flatSol.x = initSol.x;
flatSol.y = initSol.y;

flatSol.y(3,:) = 0.*flatSol.y(3,:);
flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);

Sols{1} = [flatSol.x; flatSol.y];
% gammaSols{1} = [L.*flatSol.x; ones(1, length(flatSol.x))];
foundationSols{1} = [initSol.y(1,:); y0.*ones(1, length(initSol.x))];
% cellPositions{1} = taggedPositions;

% First non-trivial solution
gammaIncremental = (gammaOld - 1)./(parameters.dt);
Sols{2} = [initSol.x; initSol.y; gammaIncremental; parameters.A.*ones(1, length(initSol.x))];
% gammaSols{2} = [L.*initSol.x; parameters.gamma];
foundationSols{2} = [parameters.Px; parameters.Py];

% cellPositions{2} = taggedPositions;

tic

% Update the solutions in time
for i = 3:numSols
    
    parameters.t = times(i);
    
    % Update the solution
    [solNew, gammaNew, EbNew, KNew, PxNew, PyNew, ANew] = UpdateAgeingAndPreSloughingSolution(solOld, parameters, solOptions);
    
    % Update the solutions, gamma, and the foundation shape
    gammaOld = parameters.gamma;
    gammaInc = (gammaNew - gammaOld)./(dt*gammaOld);
    gammaIncremental = interp1(solOld.x, gammaInc, solNew.x);
    
    parameters.gamma = interp1(solOld.x, gammaNew, solNew.x);
    %     parameters.gamma = gammaNew;
    
    parameters.Px = interp1(solOld.x, PxNew, solNew.x);
    parameters.Py = interp1(solOld.x, PyNew, solNew.x);
    
    AOld = parameters.A;
    AOld = interp1(solOld.x, AOld.*ones(1, length(solOld.x)), solNew.x);
        
    parameters.A = interp1(solOld.x, ANew.*ones(1, length(solOld)), solNew.x);
    
    parameters.currentArcLength = cumtrapz(solNew.y(1,:), parameters.gamma);
    
    sloughedAmount = trapz(solNew.y(1,:), (parameters.A > parameters.As).*(AOld.*W(parameters.currentArcLength, parameters.sigma) < 0.05));
    
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if ( (~isempty(InterX([solNew.y(2,:); solNew.y(3,:)])))||(sloughedAmount > 0) ) %
        
        Sols = Sols(1:(i - 1));
        times = times(1:(i - 1));
        foundationSols = foundationSols(1:(i - 1));
        
        break
    end
    
    
    solOld = solNew;
    
    if (length(parameters.A) > 1)
        Sols{i} = [solOld.x; solOld.y; gammaIncremental; parameters.A];
    else
        Sols{i} = [solOld.x; solOld.y; gammaIncremental; parameters.A.*ones(1, length(solOld.x))];
    end
    foundationSols{i} = [parameters.Px; parameters.Py];
    
end

toc

%% Simulate with non-zero sloughing

H = parameters.H;
T = parameters.T;
g = parameters.g;
dt = parameters.dt;
LOld = parameters.L;
gammaOld = parameters.gamma;
AOld = parameters.A;

gammaNew = gammaOld + g*dt;
parameters.gamma = gammaNew;

parameters.currentArcLength = gammaNew.*solOld.y(1,:);

sloughedAmount = trapz(times, 0.25*g*(1 + tanh(H*(times - T))));

% Define new functionx
sloughLfunct = @(l) sloughedAmount - (parameters.currentArcLength(end) - interp1(solOld.y(1,:), parameters.currentArcLength, l));

LNew = fsolve(sloughLfunct, 0, optimset('Display','off'));

parameters.L = LNew;

% Initialise the times
TMax = 7.0;
newTimes = (times(end) + dt):dt:TMax;
numSols = length(newTimes);

% Initialise new solution structures
sloughingSols = cell(numSols, 1);
gammaSloughingSols = cell(numSols, 1);
foundationSloughingSols = cell(numSols, 1);
% cellSloughingPositions = cell(numSols, 1);

sloughSolOld.x = [solOld.x, 1 2];
sloughSolOld.y = [solOld.y, repmat(solOld.y(:, end), [1 2])];

% Need to modify the solution components individually
sloughSolOld.y(1, 1:(end - 2)) = solOld.y(1, :)./(solOld.y(1, end)).*LNew;
sloughSolOld.y(1, end - 1) = LNew; % Impose new sloughing region

PxOld = parameters.Px;
parameters.Px = [PxOld, repmat(PxOld(end), [1 2])];

PyOld = parameters.Py;
parameters.Py = [PyOld, repmat(PyOld(end), [1 2])];

% parameters.A = [AOld, repmat(AOld(end), [1 2])];

parameters.currentArcLength = sloughSolOld.y(1,:).*parameters.gamma;
% parameters.gamma = [gammaOld, repmat(gammaOld(end), [1 2])];

%%
tic

% Update the solutions in time
for i = 1:numSols
    
    parameters.t = newTimes(i);
    
    % Update the solution
    [sloughSolNew, gammaNew, EbNew, KNew, PxNew, PyNew, ANew, LNew] = UpdateAgeingAndSloughingSolution(sloughSolOld, parameters, solOptions);
    
    
    % Update the solutions, gamma, and the foundation shape
    parameters.gamma = InterpolateToNewMesh(sloughSolOld.x, gammaNew, sloughSolNew.x);
    %     parameters.gamma = gammaNew;
    parameters.Px = InterpolateToNewMesh(sloughSolOld.x, PxNew, sloughSolNew.x);
    parameters.Py = InterpolateToNewMesh(sloughSolOld.x, PyNew, sloughSolNew.x);
    parameters.L = LNew;
    parameters.A = ANew;
    parameters.currentArcLength = sloughSolNew.y(1,:).*parameters.gamma;
    
    sloughSolOld = sloughSolNew;
    
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if ( (~isempty(InterX([sloughSolNew.y(2,:); sloughSolNew.y(3,:)])))) % Stopping criterion
        
        sloughingSols = sloughingSols(1:(i - 1));
        %         gammaSloughingSols = gammaSloughingSols(1:(i - 1));
        newTimes = newTimes(1:(i - 1));
        foundationSloughingSols = foundationSloughingSols(1:(i - 1));
        
        break
    end
    
    sloughingSols{i} = [sloughSolOld.x; sloughSolOld.y; gammaIncremental; parameters.A];
    foundationSloughingSols{i} = [parameters.Px; parameters.Py];
    
    %     % Find where the tagged cells are
    %     taggedPositions = zeros(length(trackedPoints),2);
    %     for j = 1:length(trackedPoints)
    %         solAtMesh = deval(sloughSolOld, solMesh(1:end - 1), [1:3]);
    %         taggedPositions(j,:) = solAtMesh(2:end, find(solAtMesh(1,:) == trackedPoints(j), 1))';
    %     end
    %
    %     cellSloughingPositions{i} = taggedPositions;
    
end

toc

%%
% sloughingSols = sloughingSols(1:(i - 1));
% gammaSloughingSols = gammaSloughingSols(1:(i - 1));
% newTimes = newTimes(1:(i - 1));
% foundationSloughingSols = foundationSloughingSols(1:(i - 1));

fullSols = [Sols; sloughingSols];
fullGammaSols = [Sols; gammaSloughingSols];
fullFoundationSols = [Sols; foundationSloughingSols];
fullTimes = [times, newTimes];
% fullCellPositions = [cellPositions; cellSloughingPositions];


%% Save the solutions
outputDirectory = '../../Solutions/AgeingAndSloughing/';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_current_sigma_w_timestepsloughing_T_2p5';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'fullSols') % Solutions
save([outputDirectory, 'gamma_', outputValues,'.mat'], 'fullGammaSols') % Gamma
save([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'fullFoundationSols') % Foundation stresses
% save([outputDirectory, 'taggedcells_', outputValues,'.mat'], 'fullCellPositions') % Migrating cells
save([outputDirectory, 'times_', outputValues, '.mat'], 'fullTimes') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

%%
outputDirectory = '../../Solutions/AgeingAndSloughing/';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_homoggrowth_timestepsloughing_T_2p5';
load([outputDirectory, 'sols_', outputValues, '.mat'], 'fullSols') % Solutions
load([outputDirectory, 'gamma_', outputValues,'.mat'], 'fullGammaSols') % Gamma
load([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'fullFoundationSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'fullTimes') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
%%
movieObj =  VideoWriter('homoggrowth_rampsloughingtime_T_2p5.avi');
movieObj.FrameRate = 15;
open(movieObj);

for i = 2:length(fullSols)
    
    figure(1)
    set(gcf, 'color', 'w');
    
    clf
    %     subplot(1, 2, 1)
    %     hold on
    %     plot(fullTimes, parameters.g.*ones(1, length(fullTimes)), 'b', 'linewidth', 2)
    %     plot(fullTimes, 0.25*(1 + tanh(parameters.H*(fullTimes - parameters.T))), 'r', 'linewidth', 2)
    %     plot(fullTimes(i).*ones(1, length(0:1e-2:1.25)), 0:1e-2:1.25, 'k', 'linewidth', 2)
    %     title('Growth vs death rate')
    %     xlabel('time t')
    %     ylim([0 1.25])
    %     xlim([0 fullTimes(end)])
    %     set(gca, 'linewidth', 1.5)
    
    
    %     subplot(1, 2, 2)
    hold on
    plot(-fullSols{i}(3,:), -fullSols{i}(4,:), 'k', 'linewidth', 2)
    plot(fullSols{i}(3,:), -fullSols{i}(4,:), 'k', 'linewidth', 2)
    %     for j = 1:3
    %         if (fullCellPositions{i}(j,1) <= 0.5)
    %             plot(fullCellPositions{i}(j, 1), -fullCellPositions{i}(j, 2), 's', 'Markersize', 15)
    %         end
    %     end
    ylim([-1.75 0.1])
    xlim([-1 1])
    title('Rod shape')
    set(gca, 'linewidth', 1.5)
    xlabel('x')
    ylabel('y')
    
    currentFrame = getframe(gcf);
    writeVideo(movieObj, currentFrame);
    
end
%
% for i = 1:length(sloughingSols)
%
%     figure(1)
%     clf
%     subplot(1, 2, 1)
%     hold on
%     plot(fullTimes, g.*ones(1, length(fullTimes)), 'b', 'linewidth', 2)
%     plot(fullTimes, 0.5*(1 + tanh(H*(fullTimes - T))), 'r', 'linewidth', 2)
%     plot(fullTimes(i).*(0:1e-2:1.5), 0:1e-2:1.5, 'k', 'linewidth', 2)
%
%     subplot(1, 2, 2)
%     hold on
%     plot(-sloughingSols{i}(3,:), -sloughingSols{i}(4,:), 'r', 'linewidth', 2)
%     plot(sloughingSols{i}(3,:), -sloughingSols{i}(4,:), 'r', 'linewidth', 2)
%     ylim([-1 0.1])
%     xlim([-1 1])
%     set(gca, 'linewidth', 1.5)
%
%     currentFrame = getframe;
%     writeVideo(movieObj, currentFrame);
%
% end

close(movieObj);
