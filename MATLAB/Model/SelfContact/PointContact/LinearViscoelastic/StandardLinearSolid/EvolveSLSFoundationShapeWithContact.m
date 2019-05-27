function EvolveSLSFoundationShapeWithContact
% Load the solutions previously computed
outputValues = 'Eb_1_beta_0p5_nu_1_kf_0p01_L0_0p125_homoggrowth';
outputDirectory = '../../../../../Solutions/LinearViscoelasticFoundation/StandardLinearSolid/';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'stresses_', outputValues,'.mat'], 'stressSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

% Initialise the current solution, gamma and the stresses

solFromData.x = Sols{end}(1,:);
solFromData.y = Sols{end}(2:end,:);

gammaFromData = 1 + (parameters.g)*times(end);
PFromData = stressSols{end}(2,:);

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

initP = PFromData([1:contIndex,contIndex:end]);

parameters.P = initP;
parameters.gamma = gammaFromData;

%%

% Solve the new contact system

% Define the ODEs and BCs
DerivFun = @(x, M, region) StandardLinearSolidFoundationContactHalfIntervalOdes(x, M, region, initSol, parameters);

% Set the boundary conditions
BcFun = @(Ml, Mr) SelfPointContactHalfIntervalBCs(Ml, Mr, parameters);

maxPoints = 1e4;

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');

tic
% Solve the system.
contactSolOld = bvp4c(DerivFun, BcFun, initSol, solOptions);

toc

%% Interpolate the new solutions on the new meshes

dt = parameters.dt;
nu = parameters.nu;
beta = parameters.beta;
K = parameters.K;
y0 = parameters.y0;

POld = parameters.P;

SOld = initSol.y(1,:);
XOld = initSol.y(2,:);
YOld = initSol.y(3,:);

DeltaOld = sqrt((XOld - SOld).^2 + YOld.^2);

SNew = deval(contactSolOld, initSol.x, 1);
XNew =  deval(contactSolOld, initSol.x, 2);
YNew =  deval(contactSolOld, initSol.x, 3);

DeltaNew = sqrt((XNew - SNew).^2 + (YNew).^2);

if (nu == 0)
    
    PNew = (1 - beta)*K.*(DeltaNew - y0);
    
else
    
    PNew = (1 - dt*beta/nu)^(-1).*(POld + dt*beta*(1 - beta)*K/nu.*(DeltaNew - y0) ...
        + (K/nu).*(DeltaNew - DeltaOld));
    
end


parameters.P = InterpolateToNewMesh(initSol.x, PNew, contactSolOld.x);


%% Run the simulation for an extra day
tMax = 2.0;
tMin = 1.0;
dt = 0.05;
parameters.dt = dt;
newTimes = tMin:dt:tMax;
numSols = length(newTimes);
contactSols = cell(numSols, 1);
stressContactSols = cell(numSols, 1);

contactSols{1} = [initSol.x; initSol.y];
stressContactSols{1} = [solFromData.y(1,:); PFromData];


contactSols{2} = [contactSolOld.x; contactSolOld.y];
stressContactSols{2} = [contactSolOld.x; parameters.P];

% momentAtContact(1) = contactSolOld.y(7, find(contactSolOld.x == 1, 1));

tic
for i = 3:numSols
    
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
    
    [contactSolNew, EbNew, gammaNew, PNew] = UpdateSLSFoundationShapeWithContact(contactSolOld, solGuess, parameters, solOptions);
    
    % Interpolate the bending stiffness to the new mesh
    %     parameters.Eb =  InterpolateToNewMesh(contactSolOld.x, EbNew, contactSolNew.x);
    
    POld = parameters.P;
    
%     if (norm(PNew - POld, 2) < 1e-4)
%         
%        contactSols = contactSosl(1:(i - 1));
%        stressContactSols = stressContactSols(1:(i - 1));
%        
%         break
%     end
    
    % Update the foundation shape
    parameters.P = InterpolateToNewMesh(contactSolOld.x, PNew, contactSolNew.x);
    
    contactSols{i} = [contactSolNew.x; contactSolNew.y];
    stressContactSols{i} = [contactSolNew.x; parameters.P];
    
    contactSolOld = contactSolNew;
    parameters.gamma = gammaNew;
    
end

toc

%%
contactSols = contactSols(1:(i - 1));
newTimes = newTimes(1:(i - 1));
% momentAtContact = momentAtContact(1:(end - 1));
stressContactSols = stressContactSols(1:(i - 1));

%% Save the solutions
outputDirectory = '../../../../../Solutions/LinearViscoelasticFoundation/StandardLinearSolid/';
outputValues = 'Eb_1_beta_0p5_nu_1_kf_0p01_L0_0p125_relaxation';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'contactSols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
save([outputDirectory, 'stresses_', outputValues,'.mat'], 'stressContactSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
% save([outputDirectory, 'contactMoments_', outputValues,'.mat'], 'momentAtContact') % Foundation stresses
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

