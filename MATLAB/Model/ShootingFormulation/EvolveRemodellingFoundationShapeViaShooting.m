function EvolveRemodellingFoundationShapeViaShooting
% Set the parameters
kf = 0.01; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
% L = 2*sqrt(3)*L0/h; %Dimensionless length
L = 1;
y0 = 0*L;
% K = kf*h/(12*w); % Dimensionless foundation stiffness
% K = 12*kf*L0^4/(w*h^3);
K = 50;
Es = 1; % Stretching stiffness
b1 = 0.5; % Bending stiffness
dt = 0.05; % Time step

% Get the initial solution from AUTO
% solData = load('../../../Data/planarmorphorodsinextkf0p01L1sol_1');
solData = load('../../../Data/seashellsK50L1sol_1');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

solFromData.y(3,:) = y0 + solFromData.y(3,:);

% sigma = w/L0; % "Width" of Wnt gradient
% w0 = 0.0*sigma;
sigma = 0.005;
w0 = 0.1;

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

solTol = 1e-4; % Set the tolerance for shooting guesses

solMesh = 0:2.5*1e-3:1;
%% Solve the initial bvp to obtain a structure for the first solution.
SOld = solFromData.y(1,:);

firstGamma = 1 + g*dt;
parameters.gamma = firstGamma;

parameters.Eb = 1 - b1.*W(SOld, sigma);

% Initialise intrinsic curvature
parameters.uHat = zeros(1, length(solFromData.x));

% Initialise foundation shape
parameters.Px = SOld;
parameters.Py = y0.*ones(1, length(solFromData.x));

% Define the ODEs and BCs
DerivFun = @(x, M) RemodellingFoundationOdes(x, M, solFromData, parameters);

% Set the tolerances and max. number of mesh points
solOptions = odeset('RelTol', 1e-4,'AbsTol', 1e-4, 'Vectorized', 'On');

tic
% Solve for the initial solution
initGuess = solFromData.y([3, 4, 7], 1)';
numSol = ShootForPreContactSolution(solFromData.y([3, 4, 7], 1)', solFromData, parameters, solOptions, solTol);

toc

initSol = GetPreContactSolution(numSol, solMesh, solFromData, parameters, solOptions);
%% Now let's iterate until the solution hits self-contact

% Update the initial data
solOld = initSol;

initS = solOld.y(1,:);
initY = solOld.y(3,:);

% Interpolate all heterogeneous quantities to the new initial mesh
parameters.currentArcLength = initS.*(parameters.gamma);

parameters.Px = initS;
parameters.Py = dt*(parameters.nu).*initY;

uHatOld = parameters.uHat;
parameters.uHat = interp1(solFromData.x, uHatOld, solMesh);

EbOld = parameters.Eb;
parameters.Eb = interp1(solFromData.x, EbOld, solMesh);

TMax = 1000.0;
times = 0:dt:TMax;
numSols = length(times);

% Initialise the solutions
Sols = cell(numSols, 1);

Sols{1} = [initSol.x; initSol.y; parameters.Px; parameters.Py; parameters.currentArcLength; parameters.Eb];

tic 
for i = 2:numSols
    
    parameters.t = times(i);
    
    % Solution guess is based on continuation of current and previous
    % solutions
    if (i == 2)
        solGuess = Sols{i - 1}([4, 5, 8], 1)' + parameters.dt*(Sols{i - 1}([4, 5, 8], 1)');
    else
        solGuess = Sols{i - 1}([4, 5, 8], 1)' + ...
            parameters.dt*(Sols{i - 1}([4, 5, 8], 1)' - Sols{i - 2}([4, 5, 8], 1)');
    end
    
    % Update the solution
    [solNew, gammaNew, EbNew, KNew, PxNew, PyNew, uHatNew] = ...
        UpdatePreContactSolutionWithShooting(solGuess, solOld, solMesh, parameters, solOptions, solTol);
    
    %     Update the incremental growth
    parameters.gamma = gammaNew;
    
    % Update the foundation shape
    parameters.Px = interp1(solOld.x, PxNew, solNew.x);
    parameters.Py = interp1(solOld.x, PyNew, solNew.x);
    
    % Update the intrinsic curvature
    parameters.uHat = interp1(solOld.x, uHatNew, solNew.x);
    
    % Update the bending stiffness
    parameters.Eb = interp1(solOld.x, EbNew, solNew.x);
    
    % Update the current arclength
    parameters.currentArcLength = cumtrapz(solNew.y(1,:), parameters.gamma.*ones(1, length(solNew.x)));
    
    % Stop the solution if the curve self-intersects
    if (HasRodHitSelfContact(solNew, parameters) )
        
        Sols = Sols(1:(i - 1));
        times = times(1:(i - 1));
        
        break
    end
    
    solOld = solNew;
    
    Sols{i} = [solOld.x; solOld.y; parameters.Px; parameters.Py; parameters.currentArcLength; parameters.Eb];
    
end

toc

%%

% Save the solutions
outputDirectory = '../../Solutions/Shooting/';
% outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_current_sigma_2w_mechanoswitchmultiplicativegrowth_linearsensitivity_mu3_0p02_h_10_T_1';
% outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_current_sigma_2w_normalised_bimodal_split_2p35w';
outputValues = 'Eb_0p5_sigmaE_0p005_nu_10_K_50_L_1_homoggrowth_w0_0p1';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% save([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
% save([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
