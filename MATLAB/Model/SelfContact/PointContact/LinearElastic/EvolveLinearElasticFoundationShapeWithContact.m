function EvolveLinearElasticFoundationShapeWithContact
% Load the solutions previously computed
outputValues = 'Eb_0p75_sigmaE_2w_k_0p02_L0_0p125_sigma_2w_area_1_mu_0_inext';
outputDirectory = '../../../../Solutions/LinearElasticFoundation/';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times

% Initialise the current solution, gamma and the stresses
solFromData.x = Sols{end}(1,:);
solFromData.y = Sols{end}(2:end,:);

gammaFromData = gammaSols{end}(2,:);

parameters.gamma = gammaFromData;

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
dt = 5*1e-2;
w0 = 0; % Width of contact region

sigma = 2*sqrt(3)*2*w/h; % "Width" of Wnt gradient
% Define the Wnt function
W = @(S, width) exp(-((S - 0.5*L)/width).^2);

eta = 1.0/trapz(solMesh, W(solMesh, sigma)); % Define eta such that the area is unit one.
mu = 0;
% 
parameters.K = K;% Foundation stiffness
parameters.L = 0.5*L; % Rod length
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

%% Construct the solution and initial guess so that bvp4c is happy

% Initialise the guesses for the contact points
possibleContactPoints = find(abs(solFromData.y(2,:) - 0.5*L) < 0.125);
contIndexOne = possibleContactPoints(1);
contIndexTwo = possibleContactPoints(end);
sC1Guess = solFromData.y(1, contIndexOne);
pCGuess = 0.0;
fCGuess = 0;

initSol.x = [solFromData.x(1:contIndexOne), solFromData.x(contIndexOne:end)];

AGuess = zeros(1, length(solFromData.x));

for i = (contIndexOne + 1):length(solFromData.x)
    AGuess(i) = trapz(solFromData.x(contIndexOne:i), ...
                    solFromData.y(2,contIndexOne:i).*(parameters.gamma(contIndexOne:i)).*sin(solFromData.y(6, contIndexOne:i)));

end


initSol.y = [solFromData.y(:, 1:contIndexOne), solFromData.y(:, contIndexOne:end);
    0.*(1:contIndexOne), AGuess(contIndexOne:end); ...
    sC1Guess.*ones(1, length(initSol.x)); ...
    pCGuess.*ones(1, length(initSol.x)); ...
    fCGuess.*ones(1, length(initSol.x))];
                
initGamma = [gammaFromData(1:contIndexOne), gammaFromData(contIndexOne:end)];
                             
initSol.y = initSol.y(:, 1:((end)/2));
initSol.y(1:2, end) = [0.5*L; 0.5*L];
initSol.y(5:6, end) = [0; 0];


initSol.x = [linspace(0, 1, contIndexOne), ...
                linspace(1, 2, length(initSol.y(1,:)) - contIndexOne)];

solMesh = initSol.x;

parameters.Eb = 1 - Eb.*W(initSol.y(1,:), sigma); % Bending stiffness`
% 
FOld = initSol.y(4,:);
GOld = initSol.y(5,:);
thetaOld = initSol.y(6,:);
n3Old = FOld.*cos(thetaOld) + GOld.*sin(thetaOld);

initGamma = initGamma(1:(end)/2);
firstGamma = initGamma.*(1 + dt*(W(initSol.y(1,:), sigma) + mu.*(n3Old - n3s)));
parameters.gamma = firstGamma;

%% Solve the initial solution

% Define the ODEs and BCs
DerivFun = @(x, M, region) LinearElasticFoundationContactOdes(x, M, region, initSol, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) SelfPointContactNormalPressureBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-5,'AbsTol', 1e-5, 'NMax', 1e6, 'Vectorized', 'On');

tic
% Solve the system. 
newSol = bvp4c(DerivFun, BcFun, initSol, solOptions);
% 
% numSol.x = solMesh;
% numSol.y = deval(numSol, solMesh);

toc


