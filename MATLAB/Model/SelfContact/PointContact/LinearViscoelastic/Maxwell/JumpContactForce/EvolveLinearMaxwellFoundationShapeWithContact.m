function EvolveLinearMaxwellFoundationShapeWithContact
% Load the solutions previously computed
outputValues = 'Eb_0p75_sigmaE_2w_nu_0p16_k_0p02_L0_0p125_sigma_2w_etaK_0p16_initforce';
outputDirectory = '../../../../../Solutions/LinearViscoelasticFoundation/Maxwell/';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'maxwellstresses_', outputValues,'.mat'], 'stressSols') % Foundation stresses
load([outputDirectory, 'intrinscurvs_', outputValues,'.mat'], 'uHatSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

% Initialise the current solution, gamma and the stresses
solFromData.x = Sols{end}(1,:);
solFromData.y = Sols{end}(2:end,:);

gammaFromData = gammaSols{end}(2,:);
PFromData = stressSols{end}(2,:);
uHatFromData = uHatSols{end}(2,:);

% parameters.gamma = initGamma(2,:);
% parameters.P = initP(2,:);

% Set the solver options for bvp4c
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

% Set the parameters
kf = 0.16; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
L = 2*sqrt(3)*L0/h; %Dimensionless length
K = kf*h/(12*w); % Dimensionless foundation stiffness
n3s = 0; % Target axial tension
Es = 1; % Stretching stiffness
Eb = 0.75; % Bending stiffness
dt = 5*1e-2;
w0 = 0; % Width of contact region

parameters.w0 = w0;
parameters.L = 0.5*L;

sigma = 2*sqrt(3)*2*(w/h); % "Width" of Wnt gradient
% Define the Wnt function
W = @(S, width) exp(-(0.5*L.*(0.5*S - 1)/width).^2);


%% Construct the solution and initial guess so that bvp4c is happy

% Initialise the guesses for the contact points
possibleContactPoints = find(abs(solFromData.y(2,:) - 0.5*L) < 0.035);
contIndexOne = possibleContactPoints(1);
sC1Guess = solFromData.y(1, contIndexOne);
fCGuess = 0;

initSol.x = [solFromData.x(1:contIndexOne), solFromData.x(contIndexOne:end)];

initSol.y = [solFromData.y(:, 1:contIndexOne), solFromData.y(:, contIndexOne:end); ...
    sC1Guess.*ones(1, length(initSol.x)); ...
    fCGuess.*ones(1, length(initSol.x))];

initGamma = [gammaFromData(1:contIndexOne), gammaFromData((contIndexOne):end)];
            
initP = [PFromData(1:contIndexOne), PFromData((contIndexOne):end)];
            
inituHat = [uHatFromData(1:contIndexOne), uHatFromData(contIndexOne:end)];

initSol.y = initSol.y(:, 1:((end - 1)/2));
initSol.y(1:2, end) = [0.5*L; 0.5*L];
initSol.y(5:6, end) = [0; 0];

initP = initP(1:((end - 1)/2));
parameters.P = initP;

inituHat = inituHat(1:((end - 1)/2));
parameters.uHat = inituHat;

initSol.x = [linspace(0, 1, contIndexOne), ...
                linspace(1, 2, length(initSol.y(1,:)) - contIndexOne)];

solMesh = initSol.x;

parameters.Eb = 1 - Eb.*W(initSol.x, sigma); % Bending stiffness`
initGamma = initGamma(1:(end - 1)/2);
firstGamma = initGamma.*(1 + (dt)*(W(initSol.x, sigma)));
parameters.gamma = firstGamma;

%% Solve the initial solution

% Define the ODEs and BCs
DerivFun = @(x, M, region) SimplifiedLinearMaxwellFoundationContactOdes(x, M, region, initSol, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) SelfPointContactBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

tic
% Solve the system. 
newSol = bvp4c(DerivFun, BcFun, initSol, solOptions);
% 
numSol.x = solMesh;
numSol.y = deval(newSol, solMesh);

toc

%%
solOld = numSol;
    
initY = numSol.y(3,:);

nu = parameters.nu;
parameters.P = initY - initSol.y(3,:) + (1 - dt*nu).*initP; 

etaK = parameters.etaK;
parameters.uHat = inituHat + etaK*dt.*initSol.y(7,:)./parameters.Eb;

% Interpolate the parameters to the new mesh
firstGamma = InterpolateToNewMesh(numSol.x, initSol.x, parameters.gamma);
firstEb = InterpolateToNewMesh(numSol.x, initSol.x, parameters.Eb);
firstP = InterpolateToNewMesh(numSol.x, initSol.x, parameters.P);
firstuHat = InterpolateToNewMesh(numSol.x, initSol.x, parameters.uHat);

parameters.gamma = firstGamma;
parameters.Eb = firstEb;
parameters.P = firstP;
parameters.uHat = firstuHat;

% Set the times we want to solve the problem for
dt = 2.5*1e-2;
parameters.dt = dt; 

TMax = 0.5;
contactTimes = 0:dt:TMax;
numSols = length(contactTimes);

% Initialise the solutions
contactSols = cell(numSols, 1);
contactGammaSols = cell(numSols, 1);
contactStressSols = cell(numSols, 1);
contactuHatSols = cell(numSols, 1);

%  The first solution is always flat

contactSols{1} = [initSol.x; initSol.y];
contactGammaSols{1} = [initSol.y(1,:); initGamma];
contactStressSols{1} = [initSol.y(1,:); initP];
contactuHatSols{1} = [initSol.y(1,:); inituHat];


% First non-trivial solution
contactSols{2} = [numSol.x; numSol.y];
contactGammaSols{2} = [numSol.y(1,:); parameters.gamma];
contactStressSols{2} = [numSol.y(1,:); parameters.P];
contactuHatSols{2} = [numSol.y(1,:), parameters.uHat];

%%
tic

% Update the solutions in time
for i = 3:numSols
        
    % Update the solution
    [solNew, gammaNew, PNew, uHatNew] = UpdateSimplifiedLinearMaxwellContactSolution(solMesh, solOld, W, parameters, solOptions);     
    
    % Update the solutions, gamma, and the spring stresses    
	parameters.gamma = InterpolateToNewMesh(solNew.x, solOld.x, gammaNew);
	parameters.P = InterpolateToNewMesh(solNew.x, solOld.x, PNew);
	parameters.uHat = InterpolateToNewMesh(solNew.x, solOld.x, uHatNew);
	parameters.Eb = InterpolateToNewMesh(solNew.x, solOld.x, parameters.Eb);
    
    figure(1)
    hold on
    plot(solNew.x, parameters.Eb)
        
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if ( (trapz(solNew.x, gammaNew) < 1)||(solNew.y(9,1) < 0) )

        contactSols = contactSols(1:(i - 1));
        contactGammaSols = contactGammaSols(1:(i - 1));
        contactTimes = contactTimes(1:(i - 1));
        contactStressSols = contactStressSols(1:(i - 1));
        contactuHatSols = contactuHatSols(1:(i - 1));
        
        break
    end
    
    solOld = solNew;
    
    contactSols{i} = [solOld.x; solOld.y];
    contactGammaSols{i} = [solOld.y(1,:); parameters.gamma];
    contactStressSols{i} = [solOld.y(1,:); parameters.P];
    contactuHatSols{i} = [solOld.y(1,:); parameters.uHat];
                        
end

toc

% Save the solutions


outputValues = 'postcontact_Eb_0p75_sigmaE_2w_nu_0p16_k_0p02_L0_0p125_sigma_2w_etaK_0p16_initforce';
outputDirectory = '../../../../../Solutions/LinearViscoelasticFoundation/Maxwell/';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'contactSols') % Solutions
save([outputDirectory, 'gamma_', outputValues,'.mat'], 'contactGammaSols') % Gamma
save([outputDirectory, 'maxwellstresses_', outputValues,'.mat'], 'contactStressSols') % Foundation stresses
save([outputDirectory, 'intrinscurvs_', outputValues,'.mat'], 'contactuHatSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'contactTimes') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
