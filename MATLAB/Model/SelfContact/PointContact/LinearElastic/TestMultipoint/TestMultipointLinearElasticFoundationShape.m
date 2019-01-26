function TestMultipointLinearElasticFoundationShape
% Load the solutions previously computed
outputValues = 'Eb_0p75_sigmaE_2w_k_0p02_L0_0p125_sigma_2w_area_1_mu_0_inext';
outputDirectory = '../../../../../Solutions/LinearElasticFoundation/';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times

% Initialise the current solution, gamma and the stresses
solFromData.x = Sols{end - 1}(1,:);
solFromData.y = Sols{end - 1}(2:end,:);

gammaFromData = gammaSols{end - 1}(2,:);

solMesh = solFromData.x;

% Set the solver options for bvp4c
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

% Set the parameters
kf = 0.16; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
L = 2*sqrt(3)*L0/h; %Dimensionless length
y0 = 0;
K = kf*h/(12*w); % Dimensionless foundation stiffness
n3s = 0; % Target axial tension
Es = 1; % Stretching stiffness
Eb = 0.75; % Bending stiffness
dt = 2.5*1e-2;
w0 = 0; % Width of contact region

sigma = 2*sqrt(3)*2*w/h; % "Width" of Wnt gradient
% Define the Wnt function
W = @(S, width) exp(-(L*(S - 0.5)/width).^2);

eta = 1.0/trapz(solMesh, W(solMesh, sigma)); % Define eta such that the area is unit one.
mu = 0;
% 
parameters.K = K;% Foundation stiffness
parameters.L = L; % Rod length
parameters.y0 = y0;
parameters.sigma = sigma; % Width of wnt gradient
parameters.mu = mu; % Rate of mechanical inhibition
parameters.eta = eta; % Rate of chemical change
parameters.n3s = n3s; % Target axial stress
parameters.Es = Es; % Stretch stiffness
% parameters.Eb = 1 - Eb.*W(solMesh, sigma); % Bending stiffness
parameters.ext = 0; % Exstensibility
parameters.dt = dt; % Time step
parameters.w0 = w0;

%% Initialise a guess for the contact solution

% We need to duplicate the midpoint so that bvp4c can read it in correctly
initSol.x = [solFromData.x(1:300), 0.5, 0.5, solFromData.x(301:end)];
initSol.y = [solFromData.y(:, 1:300), ...
                    [0.5*L; 0.5*L; 0.5*(solFromData.y(3,300) + solFromData.y(3,301)); ...
                    0.5*(solFromData.y(4,300) + solFromData.y(4,301)); 0; 0; 
                     0.5*(solFromData.y(7,300) + solFromData.y(7,301))], ...
                     [0.5*L; 0.5*L; 0.5*(solFromData.y(3,300) + solFromData.y(3,301)); ...
                    0.5*(solFromData.y(4,300) + solFromData.y(4,301)); 0; 0; 
                     0.5*(solFromData.y(7,300) + solFromData.y(7,301))], ...
                    solFromData.y(:, 301:end)];
                
solMesh = initSol.x;

FOld = initSol.y(4,:);
GOld = initSol.y(5,:);
thetaOld = initSol.y(6,:);
n3Old = FOld.*cos(thetaOld) + GOld.*sin(thetaOld);

initGamma = [gammaFromData(1:300), 0.5*(gammaFromData(300) + gammaFromData(301)), ...
                0.5*(gammaFromData(300) + gammaFromData(301)), gammaFromData(301:end)];
firstGamma = initGamma.*(1 + dt*(W(solMesh, sigma) + mu.*(n3Old - n3s)));
parameters.gamma = firstGamma;

parameters.Eb = 1 - Eb.*W(solMesh, sigma); % Bending stiffness


% Define the ODEs and BCs
DerivFun = @(x, M, region) HalfIntervalOdes(x, M, region, initSol, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) HalfIntervalBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

tic
% Solve the system. 
numSol = bvp4c(DerivFun, BcFun, initSol, solOptions);

toc
