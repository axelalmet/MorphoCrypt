function EvolveShapeWithRemodellingFoundation
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
b1 = 0; % Bending stiffness
dt = 0.01; % Time step

% Get the initial solution from AUTO
solData = load('../../../Data/planarmorphorodsinextkf0p01L1sol_1');
% solData = load('../../../Data/seashellsK50L1sol_1');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

solFromData.y(3,:) = y0 + solFromData.y(3,:);

% % Load previous heterogeneous solution to obtain a better guess for the solution
% outputDirectory = '../../Solutions/RemodellingFoundation/';
% outputValues = 'Eb_0p5_init_sigmaE_3w_nu_10_kf_0p01_L0_0p125_homoggrowth';
% load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% solFromData.x = Sols{2}(1,:);
% solFromData.y = Sols{2}(2:(end - 1), :);

sigma = w/L0; % "Width" of Wnt gradient
w0 = 0.0*sigma;
% sigma = 0.005;

% % Define the Wnt function
W = @(S, width) exp(-((S)./width).^2);

parameters.W = W;

eta = 1.0/24; % Growth timescale (24 hours)
etaF = eta^(-1);

mu3 = 0.02;
    
g = 1;
N = 2;
Q = 0;

parameters.mu3 = mu3;
parameters.g = g;
parameters.currentArcLength = solFromData.y(1,:);
parameters.K = K;% Foundation stiffness
parameters.L = 0.5*L; % Rod length
parameters.y0 = y0; % Rod centreline
parameters.sigma = sigma; % Width of wnt gradient
parameters.w0 = w0;
parameters.eta = eta; % Rate of chemical change
parameters.Es = Es; % Stretch stiffness
parameters.b1 = b1;
parameters.Eb = 1; % Bending stiffness
parameters.ext = 0; % Exstensibility
parameters.nu = 10*etaF*eta; % Foundation relaxation timescale
parameters.chi = 10*parameters.nu; % Curvature relaxation timescale
parameters.etan3 = 0*eta*etaF; % Relaxation of target stress
parameters.etaK = 0*eta*etaF; % Curvature relaxation timescale
parameters.dt = dt; % Time step
parameters.N = N;
parameters.Q = Q;

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


%% Solve the initial bvp to obtain a structure for the first solution.
SOld = solFromData.y(1,:);
FOld = solFromData.y(4,:);
GOld = solFromData.y(5,:);
thetaOld = solFromData.y(6,:);
n3Old = FOld.*cos(thetaOld) + GOld.*sin(thetaOld);

% Initialise non-unit growth
firstGamma = 1 + g*dt;
parameters.gamma = firstGamma;

% parameters.Eb = 1 - b1.*exp(-((solFromData.y(1,:))./sigma).^2);

parameters.Eb = 1;

% Initialise intrinsic curvature
parameters.uHat = zeros(1, length(solFromData.x));

% Initialise foundation shape
parameters.Px = SOld;
parameters.Py = y0.*ones(1, length(solFromData.x));

% Define the ODEs and BCs
DerivFun = @(x, M) RemodellingFoundationOdes(x, M, solFromData, parameters);
% DerivFun = @(x, M) RemodellingFoundationWithHorizontalRepulsionOdes(x, M, solFromData, parameters);

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

%%
% Update the initial data
solOld = initSol;

initS = solOld.y(1,:);
initY = solOld.y(3,:);

% Interpolate all heterogeneous quantities to the new initial mesh
% parameters.Eb = 1 - b1.*W(initS, sigma);
% firstGamma = interp1(solFromData.y(1,:), firstGamma, initSol.y(1,:));
% parameters.gamma = firstGamma;
% parameters.currentArcLength = cumtrapz(initS, parameters.gamma.*ones(1, length(initS)));
parameters.currentArcLength = initS.*(parameters.gamma);

parameters.Px = initS;
parameters.Py = dt*(parameters.nu).*initY;

uHatOld = parameters.uHat;
parameters.uHat = interp1(solFromData.x, uHatOld, initSol.x);

TMax = 1000.0;
times = 0:dt:TMax;
numSols = length(times);

% Initialise the solutions
Sols = cell(numSols, 1);
foundationSols = cell(numSols, 1);

%  The first solution is always flat
flatSol.x = initSol.x;
flatSol.y = initSol.y;

flatSol.y(3,:) = 0.*flatSol.y(3,:);
flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);

Sols{1} = [flatSol.x; flatSol.y; zeros(1, length(flatSol.x))];
foundationSols{1} = [initSol.y(1,:); y0.*ones(1, length(initSol.x))];

% First non-trivial solution
% gammaIncremental = (firstGamma - 1)/dt;
% Sols{2} = [initSol.x; initSol.y; gammaIncremental];
Sols{2} = [initSol.x; initSol.y; parameters.uHat];
% Sols{2} = [initSol.x; initSol.y; parameters.Eb];
foundationSols{2} = [parameters.Px; parameters.Py];

% Update the solutions in time

tic 
for i = 3:numSols
    
    parameters.t = times(i);
    
    solCurrent.x = Sols{i - 1}(1,:);
    solCurrent.y = Sols{i - 1}(2:end,:);
    
    solPrevious.x = Sols{i - 2}(1,:);
    solPrevious.y = Sols{i - 2}(2:end,:);
    
    solGuess = ConstructNewGuess(solCurrent, solPrevious);
    
    % Update the solution
    [solNew, gammaNew, EbNew, KNew, PxNew, PyNew, uHatNew] = UpdateRemodellingFoundationSolution(solOld, solGuess, parameters, solOptions);
    
    %     Update the incremental growth
%     gammaOld = parameters.gamma;
%     gammaInc = (gammaNew - gammaOld)./(dt*gammaOld);
%     gammaIncremental = interp1(solOld.x, gammaInc, solNew.x);
    
    % Update gamma to the new mesh
%     parameters.gamma = interp1(solOld.x, gammaNew, solNew.x);
        parameters.gamma = gammaNew;
    
    % Update the foundation shape
    parameters.Px = interp1(solOld.x, PxNew, solNew.x);
    parameters.Py = interp1(solOld.x, PyNew, solNew.x);
    
    % Update the intrinsic curvature
    parameters.uHat = interp1(solOld.x, uHatNew, solNew.x);
    
    % Update the bending stiffness
%         parameters.Eb = interp1(solOld.x, EbNew, solNew.x);
    
    % Update the current arclength
%     parameters.currentArcLength = cumtrapz(solNew.y(1,:), parameters.gamma.*ones(1, length(solNew.x)));
    
    % Stop the solution if the curve self-intersects
    if (HasRodHitSelfContact(solNew, parameters)||(times(i) > 3.0) )
        
        Sols = Sols(1:(i - 1));
        times = times(1:(i - 1));
        foundationSols = foundationSols(1:(i - 1));
        
        break
    end
    
    solOld = solNew;
    
%         Sols{i} = [solOld.x; solOld.y; gammaIncremental];
%     Sols{i} = [solOld.x; solOld.y; parameters.Eb];
    Sols{i} = [solOld.x; solOld.y; parameters.uHat];
    foundationSols{i} = [parameters.Px; parameters.Py];
    
end

toc

%%
Sols = Sols(1:(i - 1));
times = times(1:(i - 1));
foundationSols = foundationSols(1:(i - 1));

% Save the solutions
outputDirectory = '../../Solutions/RemodellingFoundation/';
% outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_current_sigma_2w_mechanoswitchmultiplicativegrowth_linearsensitivity_mu3_0p02_h_10_T_1';
% outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_current_sigma_2w_normalised_bimodal_split_2p35w';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_chi_100_expgrowth';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% save([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
save([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
