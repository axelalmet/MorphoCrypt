function EvolveShapeWithLinearElasticFoundation
% close all

% Set the parameters
kf = 0.0; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod`
% L = 2*sqrt(3)*L0/h; %Dimensionless length
L = 1;
y0 = 0;
% K = kf*h/(12*w); % Dimensionless foundation stiffness
K = 12*kf*L0^4/(w*h^3);
n3s = 0; % Target axial tension
Es = 1; % Stretching stiffness
b1 = 0;
Eb = 1; % Bending stiffness
dt = 5*1e-2;

% Get the initial solution from AUTO
solData = load('../../../Data/planarmorphorodsinextkf0L1sol_1');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

sigma = 2*sqrt(3)*w/h; % "Width" of Wnt gradient
% Define the Wnt function
W = @(S, width) exp(-((S - 0.5*S(end))./width).^2);

eta = 1/24;
mu = 0;
% 
parameters.K = K;% Foundation stiffness
parameters.L = L; % Rod length
parameters.y0 = y0;
parameters.sigma = sigma; % Width of wnt gradient
parameters.mu = mu; % Rate of mechanical inhibition
parameters.eta = eta; % Rate of chemical change
parameters.Es = Es; % Stretch stiffness
% parameters.Eb = 1 - Eb.*W(solFromData.x, sigma); % Bending stiffness
parameters.Eb = 1;
parameters.etaK = 0;
parameters.ext = 0; % Extensibility
parameters.dt = dt; % Time step

%% Solve the initial bvp to obtain a structure for the first solution.

gammaOld = 1;
% firstGamma = gammaOld.*(1 + dt*(W(solFromData.x, sigma) + mu.*(n3Old - n3s)));
firstGamma = gammaOld + dt;
parameters.gamma = firstGamma;

% Initialise intrinsic curvature
parameters.uHat = zeros(1, length(solFromData.x));

% Define the ODEs and BCs
DerivFun = @(x, M) LinearElasticFoundationOdes(x, M, solFromData, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters);
% BcFun = @(Ml, Mr) HalfIntervalBCs(Ml, Mr, parameters);

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

TMax = 50.0;
times = 0:dt:TMax;
numSols = length(times);


uHatOld = parameters.uHat;
parameters.uHat = interp1(solFromData.x, uHatOld, initSol.x);

% Initialise the solutions
Sols = cell(numSols, 1);
% gammaSols = cell(numSols, 1);
uHatSols = cell(numSols, 1);

%  The first solution is always flat
flatSol.x = initSol.x;
flatSol.y = initSol.y;  

flatSol.y(3,:) = 0.*flatSol.y(3,:);
flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);

Sols{1} = [flatSol.x; flatSol.y];
% gammaSols{1} = [L.*flatSol.x; ones(1, length(flatSol.x))];
uHatSols{1} = [flatSol.x; 0.*flatSol.x];

% First non-trivial solution
Sols{2} = [initSol.x; initSol.y];
% gammaSols{2} = [L.*initSol.x; gammaOld];
uHatSols{2} = [L.*initSol.x, parameters.uHat];

tic
% Update the solutions in time
for i = 3:numSols
    
    % Update the solution
    [solMeshNew, solNew, gammaNew, uHatNew] = UpdateLinearElasticSolution(solMesh, solOld, W, parameters, solOptions);
    
    % Update the solutions and gamma
    %     gammaOld = interp1(solOld.x, gammaNew, solNew.x);
    %     parameters.gamma = gammaOld;
    parameters.gamma = gammaNew;
    parameters.uHat = interp1(solOld.x, uHatNew, solNew.x);
    
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    %     if ( (trapz(solOld.x, gammaOld) < 1)||(~isempty(InterX([solNew.y(2,:); solNew.y(3,:)]))) )
    if (~isempty(InterX([solNew.y(2,:); solNew.y(3,:)]))) 
        
        Sols = Sols(1:(i - 1));
        %         gammaSols = gammaSols(1:(i - 1));
        times = times(1:(i - 1));
        uHatSols = uHatSols(1:(i - 1));
        
        break
    end
    
    solOld = solNew;
    solMesh = solMeshNew;
    
    Sols{i} = [solOld.x; solOld.y];
    %     gammaSols{i} = [L.*solNew.x; gammaOld];
    uHatSols{i} = [L.*solNew.x; parameters.uHat];
    
end
toc

% Save the files.   
% Sols = Sols(1:(end - 1));
% gammaSols = gammaSols(1:(end - 1));
% times = times(1:(end - 1));

outputDirectory = '../../Solutions/LinearElasticFoundation/'; 
outputValues = 'Eb_1_kf_0_L0_0p125_homoggrowth';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% save([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
save([outputDirectory, 'intrinscurvs_', outputValues,'.mat'], 'uHatSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
