function EvolveShapeWithLinearKelvinVoigtFoundation
% Set the parameters
kf = 0.01; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
% L = 2*sqrt(3)*L0/h; %Dimensionless length
L = 1;
y0  = 6*L;
% K = kf*h/(12*w); % Dimensionless foundation stiffness
K = 12*kf*L0^4/(w*h^3);
n3s = 0; % Target axial tension
Es = 1; % Stretching stiffness
b1 = 0;
dt = 2.5*1e-2; % Time step

% Get the initial solution from AUTO 
solData = load('../../../../Data/planarmorphorodsinextkf0p01L1sol_1');
% Cartesian basis

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

% Translate solution up
solFromData.y(3,:) = y0 + solFromData.y(3,:);

% 
sigma = 2*sqrt(3)*2*w/h; % "Width" of Wnt gradient
% % Define the Wnt function
W = @(S, width) exp(-(L*(S - 0.5)/width).^2);

eta = 1/24;
etaF = eta^(-1);
etaV = eta^(-1);

parameters.g = 1;
parameters.K = K;% Foundation stiffness
parameters.L = L; % Rod length
parameters.y0 = y0; % Rod centreline
parameters.sigma = sigma; % Width of wnt gradient
parameters.eta = eta; % Rate of chemical change
parameters.Es = Es; % Stretch stiffness
parameters.Eb = 1; % Bending stiffness
parameters.ext = 0; % Exstensibility
parameters.nu = 1*etaV*eta/kf; % Foundation relaxation timescale
parameters.etaK = 0*kf/(eta*etaF); % Curvature relaxation timescale
parameters.dt = dt; % Time step

%% Solve the initial bvp to obtain a structure for the first solution.

% Cartesian basis solution extraction
SOld = solFromData.y(1,:);
xOld = solFromData.y(2,:);
yOld = solFromData.y(3,:);

DeltaOld = ((xOld - SOld).^2 + yOld.^2).^0.5;

gammaOld = 1;
firstGamma = gammaOld + (parameters.g)*dt;
parameters.gamma = firstGamma;

% parameters.P = (parameters.nu)*(parameters.dt)^(-1).*(DeltaOld - (parameters.y0));
% parameters.P = yOld + (parameters.nu)*(parameters.dt)^(-1).*(yOld - (parameters.y0));
parameters.P = zeros(1, length(solFromData.x));

% Initialise intrinsic curvature
parameters.uHat = zeros(1, length(solFromData.x));

% Define the ODEs and BCs
DerivFun = @(x, M) LinearKelvinVoigtFoundationOdes(x, M, solFromData, parameters);
% DerivFun = @(x, M) SimplifiedLinearKelvinVoigtFoundationOdes(x, M, solFromData, parameters);

% Set the boundary conditions 
BcFun = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters);

tic
% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

% Solve the system. 
numSol = bvp4c(DerivFun, BcFun, solFromData, solOptions);

toc
% plot(initSol.x, initSol.y(3,:))

initSol.x = solFromData.x;
initSol.y = deval(numSol, solFromData.x); 
% initSol = numSol;
                            
%%  
% gammaOld = interp1(solFromData.x, firstGamma, initSol.x);
% parameters.gamma = gammaOld;

solOld = initSol;
initS = initSol.y(1,:);
initX = initSol.y(2,:);
initY = initSol.y(3,:);
initDelta = sqrt((initX - initS).^2 + initY.^2);

% parameters.P = parameters.nu*(dt)^(-1).*(initDelta - y0); 
parameters.P = initDelta + (parameters.nu/dt).*(initDelta - y0);

inituHat = parameters.uHat;
parameters.uHat = interp1(solFromData.x, inituHat, initSol.x);


solMesh = initSol.x;
    
% Set the times we want to solve the problem for

TMax = 50.0;
times = 0:dt:TMax;
numSols = length(times);

% Initialise the solutions
Sols = cell(numSols, 1);
% gammaSols = cell(numSols, 1);
stressSols = cell(numSols, 1);

%  The first solution is always flat
flatSol.x = initSol.x;
flatSol.y = initSol.y;  

flatSol.y(3,:) = 0.*flatSol.y(3,:) + y0;
flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);

Sols{1} = [flatSol.x; flatSol.y];
% gammaSols{1} = [L.*flatSol.x; ones(1, length(flatSol.x))];
stressSols{1} = [flatSol.x; 0.*flatSol.x];

% First non-trivial solution
Sols{2} = [initSol.x; initSol.y];
% gammaSols{2} = [L.*initSol.x; gammaOld];
stressSols{2} = [L.*initSol.x; parameters.P];


tic

% Update the solutions in time
for i = 3:numSols
        
    % Update the solution
    [solMeshNew, solNew, gammaNew, PNew, uHatNew] = UpdateLinearKelvinVoigtSolution(solMesh, solOld, W, parameters, solOptions);     
%         [solMeshNew, solNew, gammaNew, PNew, uHatNew] = UpdateSimplifiedLinearKelvinVoigtSolution(solMesh, solOld, W, parameters, solOptions);

    % Update the solutions, gamma, and the spring stresses
%     gammaOld = interp1(solOld.x, gammaNew, solNew.x);
    gammaOld = gammaNew;
    parameters.gamma = gammaOld;
    parameters.P = interp1(solOld.x, PNew, solNew.x);
%     parameters.P = PNew;
    parameters.uHat = interp1(solOld.x, uHatNew, solNew.x);
    
    
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    %     if ( (trapz(solOld.x, gammaOld) < 1)||(~isempty(InterX([solNew.y(2,:); solNew.y(3,:)]))) )
    if (~isempty(InterX([solNew.y(2,:); solNew.y(3,:)]))) 
        
        Sols = Sols(1:(i - 1));
        %         gammaSols = gammaSols(1:(i - 1));
        times = times(1:(i - 1));
        stressSols = stressSols(1:(i - 1));
        
        break
    end
    
    solOld = solNew;    

    Sols{i} = [solOld.x; solOld.y];
%     gammaSols{i} = [L.*solNew.x; gammaOld];
    stressSols{i} = [L.*solNew.x; parameters.P];
    
    solMesh = solMeshNew;
                        
end

toc

% Save the solutions

outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/KelvinVoigt/'; 
outputValues = 'Eb_1_nu_100_kf_0p01_L0_0p125_homoggrowth';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% save([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
save([outputDirectory, 'kvstresses_', outputValues,'.mat'], 'stressSols') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
% save([outputDirectory, 'intrinscurvs_', outputValues,'.mat'], 'uHatSols') % Foundation stresses
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

