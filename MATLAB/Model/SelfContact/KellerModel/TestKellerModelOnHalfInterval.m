function TestKellerModelOnHalfInterval
%% Set-up
P = 3.25;
L = pi;

parameters.P = P;
parameters.L = L;

% First derive the buckled solution
solMesh = 0:1e-2:1;

initSol.x = solMesh;
initSol.y = [L.*solMesh;
             sin(L.*solMesh + 0.5*sin(2.*L.*solMesh)); ...
             -cos(L.*solMesh + 0.5*sin(2.*L.*solMesh)); ...
             0.*solMesh - 0.5*2*cos(2.*L.*solMesh); ...
             0.*solMesh - 0.5*4*sin(2.*L.*solMesh); ...
             L.*solMesh + 0.5*sin(2.*L.*solMesh); ...
             ones(1, length(solMesh)) + 0.5*2*cos(L.*2.*solMesh)];

% Define the ODEs and BCs
DerivFun = @(x, M) KellerModelPreContactOdes(x, M, parameters);
    
% Set the boundary conditions 
BcFun = @(Ml, Mr) KellerModelHalfIntervalPreContactBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

tic
% Solve the system. 
newSol = bvp4c(DerivFun, BcFun, initSol, solOptions);

numSol.x = solMesh;
numSol.y = deval(newSol, solMesh);

toc

% figure(1)
% hold on
% plot(numSol.y(2,:), numSol.y(3,:), 'k', 'linewidth', 1.5)
% plot(-numSol.y(2,:), numSol.y(3,:), 'k', 'linewidth', 1.5)
% title('p = 3.25')
% xlim([-1 1])

%%

P = 4.75;

parameters.P = P;

% Define the ODEs and BCs
DerivFun = @(x, M) KellerModelPreContactOdes(x, M, parameters);
    
% Set the boundary conditions 
BcFun = @(Ml, Mr) KellerModelPreContactBCs(Ml, Mr, parameters);

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

tic
% Solve the system. 
newSol = bvp4c(DerivFun, BcFun, numSol, solOptions);

numSol.x = solMesh;
numSol.y = deval(newSol, solMesh);

toc

% figure(2)
% hold on
% plot(numSol.y(2,:), numSol.y(3,:), 'k', 'linewidth', 1.5)
% plot(-numSol.y(2,:), numSol.y(3,:), 'k', 'linewidth', 1.5)
% title('p = 4.75')
% xlim([-1 1])

%% Do a simple continuation step to contact

tic
parameters.P = 5.247;
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

solOld.y = [numSol.y(:, 1:(end - 1)/2), numSol.y(:, ((end + 1)/2):end)];
solOld.y = [solOld.y; 0.5*pi.*ones(1, length(solOld.y(1,:))); 0.*ones(1, length(solOld.y(1,:)))];
solOld.x = [linspace(0, 1, 0.5*length(solOld.y(1,:))), linspace(1, 2, 0.5*length(solOld.y(1,:)) + 1)];
solMesh = solOld.x;

% Define the ODEs and BCs
DerivFun = @(x, M, region) KellerModelHalfIntervalContactOdes(x, M, region, parameters);

% Set the boundary conditions
BcFun = @(Ml, Mr) KellerModelHalfIntervalContactBCs(Ml, Mr, parameters);

contactSol = bvp4c(DerivFun, BcFun, solOld, solOptions);

% contactSol.x = solMesh;
% contactSol.y = deval(contactSol, solMesh);

toc

% figure(3)
% hold on
% plot(contactSol.y(2,:), contactSol.y(3,:), 'k', 'linewidth', 1.5)
% plot(-contactSol.y(2,:), contactSol.y(3,:), 'k', 'linewidth', 1.5)
% title('p = 5.247')
% xlim([-1 1])

%% Continue the contact solution until we get close to the region contact case

contactSolOld = contactSol;

dP = 0.2;
PMax = 10.2;

tic
while (parameters.P < PMax)
    parameters.P = parameters.P + dP;
    
    % Define the ODEs and BCs
    DerivFun = @(x, M, region) KellerModelHalfIntervalContactOdes(x, M, region, parameters);
    
    % Set the boundary conditions
    BcFun = @(Ml, Mr) KellerModelHalfIntervalContactBCs(Ml, Mr, parameters);
    
    contactSolNew = bvp4c(DerivFun, BcFun, contactSolOld, solOptions);
    
%     contactSolNew.x = solMesh;
%     contactSolNew.y = deval(NewSol, solMesh);
    
    contactSolOld = contactSolNew;
%     
%     figure(5)
%     hold on
%     plot(contactSolOld.y(1,:), contactSolOld.y(4,:))
%     
%     figure(4)
%     hold on
%     plot(contactSolOld.y(1,:), contactSolOld.y(7,:))
    
    
end
toc

%%

tic
parameters.P = 10.337;
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');

% Define the ODEs and BCs
DerivFun = @(x, M, region) KellerModelHalfIntervalContactOdes(x, M, region, parameters);

% Set the boundary conditions
BcFun = @(Ml, Mr) KellerModelHalfIntervalContactBCs(Ml, Mr, parameters);

contactSol = bvp4c(DerivFun, BcFun, contactSolOld, solOptions);


toc`

figure(5)
hold on
plot(contactSol.y(2,:), contactSol.y(3,:), 'k', 'linewidth', 1.5)
plot(-contactSol.y(2,:), contactSol.y(3,:), 'k', 'linewidth', 1.5)
title('p = 10.337')
xlim([-1 1])


%% Test for contact along a region

% Construct the new solution guess

solOld = contactSol;
solOld.y(1,37:end) = solOld.y(1, 37:end) + 2*0.04;
solOld.y(3,37:end) = solOld.y(3, 37:end) + 2*0.04;
solOld.y(8,:) = 0.04;

solOld.y = solOld.y(1:8,:);

tic
parameters.P = 10.4;
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', 1e6, 'Vectorized', 'On');


% Define the ODEs and BCs
DerivFun = @(x, M, region) KellerModelHalfIntervalContactRegionOdes(x, M, region, parameters);

% Set the boundary conditions
BcFun = @(Ml, Mr) KellerModelHalfIntervalContactRegionBCs(Ml, Mr, parameters);

contactRegionSol = bvp4c(DerivFun, BcFun, solOld, solOptions);

% contactRegionSol.x = solMesh;
% contactRegionSol.y = deval(contactRegionSol, solMesh);

toc
%% Run a continuation in P to make sure it's doing what it's supposed to!

PMax = 20.0;
dP = 0.1;

contactRegionSolOld = contactRegionSol;

tic
while (parameters.P < PMax)
    parameters.P = parameters.P + dP;
    
    % Define the ODEs and BCs
    DerivFun = @(x, M, region) KellerModelHalfIntervalContactRegionOdes(x, M, region, parameters);
    
    % Set the boundary conditions
    BcFun = @(Ml, Mr) KellerModelHalfIntervalContactRegionBCs(Ml, Mr, parameters);
    
    NewSol = bvp4c(DerivFun, BcFun, contactRegionSolOld, solOptions);
%     
%     contactRegionSolNew.x = solMesh;
%     contactRegionSolNew.y = deval(NewSol, solMesh);
    
    contactRegionSolOld = NewSol;
    
    figure(1)
    hold on
    plot(contactRegionSolOld.y(1,:), contactRegionSolOld.y(3,:))

    figure(2)
    hold on
    plot(contactRegionSolOld.y(1,:), contactRegionSolOld.y(5,:))
    
    figure(3)
    hold on
    plot(contactRegionSolOld.y(1,:), contactRegionSolOld.y(4,:))
    
     
end
toc


figure(6)
hold on
plot(contactRegionSolOld.y(2,:), contactRegionSolOld.y(3,:))
plot(-contactRegionSolOld.y(2,:), contactRegionSolOld.y(3,:))
title('p = 20')
xlim([-1 1])
