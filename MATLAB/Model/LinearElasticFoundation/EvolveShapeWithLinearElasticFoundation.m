function EvolveShapeWithLinearElasticFoundation
% close all

% Set the parameters
kf = 0.01; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod`
% L = 2*sqrt(3)*L0/h; %Dimensionless length
L = 1;
y0 = 0;
% K = kf*h/(12*w); % Dimensionless foundation stiffness
K = 12*kf*L0^4/(w*h^3);
Es = 1; % Stretching stiffness
b1 = 0;
Eb = 1; % Bending stiffness
dt = 1e-2;

% Get the initial solution from AUTO
solData = load('../../../Data/planarmorphorodsinextkf0p01L1sol_1');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';


sigma = w/L0; % "Width" of Wnt gradient
w0 = 0.0*sigma;

% Define the Wnt function
W = @(S, width) exp(-((S - 0.5*S(end))./width).^2);

parameters.W = W;

eta = 1.0/24; % Growth timescale (24 hours)
etaF = eta^(-1);
etaK = etaF;

g = 1;

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
parameters.Eb = Eb; % Bending stiffness
parameters.ext = 0; % Exstensibility
parameters.chi = 0.1*eta*etaF; % Curvature relaxation timescale
parameters.dt = dt; % Time step

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

gammaOld = 1;
% firstGamma = gammaOld.*(1 + dt*(W(solFromData.x, sigma) + mu.*(n3Old - n3s)));
firstGamma = gammaOld + g*dt;
parameters.gamma = firstGamma;

% Initialise intrinsic curvature
parameters.uHat = zeros(1, length(solFromData.x));

% Define the ODEs and BCs
DerivFun = @(x, M) LinearElasticFoundationOdes(x, M, solFromData, parameters);

% Set the boundary conditions 
% BcFun = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters);
BcFun = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 2*length(solFromData.x), 'Vectorized', 'On');

tic
% Solve the system. 
numSol = bvp4c(DerivFun, BcFun, solFromData, solOptions);
toc

% plot(initSol.x, initSol.y(3,:))   

initSol = numSol;
% 
% initSol.x = solFromData.x;
% initSol.y = deval(numSol, solFromData.x); 
                            
%%  
% gammaOld = interp1(solFromData.x, firstGamma, initSol.x);
parameters.gamma = gammaOld;
solOld = initSol;
    
solMesh = numSol.x;

% Set the times we want to solve the problem for

TMax = 5.0;
times = 0:dt:TMax;
numSols = length(times);


uHatOld = parameters.uHat;
parameters.uHat = interp1(solFromData.x, uHatOld, initSol.x);

% Initialise the solutions
Sols = cell(numSols, 1);

%  The first solution is always flat
flatSol.x = initSol.x;
flatSol.y = initSol.y;  

flatSol.y(3,:) = 0.*flatSol.y(3,:);
flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);

Sols{1} = [flatSol.x; flatSol.y; zeros(1, length(flatSol.x))];

% First non-trivial solution
Sols{2} = [initSol.x; initSol.y; parameters.uHat];

tic
% Update the solutions in time
for i = 3:numSols
    
    % Update the solution
    [solNew, gammaNew, uHatNew] = UpdateLinearElasticSolution(solOld, parameters, solOptions);
    
    % Update the solutions and gamma
    %     gammaOld = interp1(solOld.x, gammaNew, solNew.x);
    %     parameters.gamma = gammaOld;
    parameters.gamma = gammaNew;
    parameters.uHat = interp1(solOld.x, uHatNew, solNew.x);
    
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if (HasRodHitSelfContact(solNew, parameters) )
        
        Sols = Sols(1:(i - 1));
        times = times(1:(i - 1));
        
        break
    end
    
    solOld = solNew;    
    Sols{i} = [solOld.x; solOld.y; parameters.uHat];
    %     gammaSols{i} = [L.*solNew.x; gammaOld];
    
end
toc

%%
% % Save the files.   
% Sols = Sols(1:);
% times = times(1:42);

outputDirectory = '../../Solutions/LinearElasticFoundation/'; 
outputValues = 'Eb_1_kf_0p01_L0_0p125_homoggrowth_chi_0p1';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% save([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
% save([outputDirectory, 'intrinscurvs_', outputValues,'.mat'], 'uHatSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
