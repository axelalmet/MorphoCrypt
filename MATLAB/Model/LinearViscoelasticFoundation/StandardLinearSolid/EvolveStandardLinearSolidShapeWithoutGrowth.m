function EvolveStandardLinearSolidShapeWithoutGrowth
% Load the solutions previously computed
outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/StandardLinearSolid/';
nuValueNames = {'0p05', '0p1', '0p2', '0p25', '0p3', '0p5', '1'};
timeSample = 1.0;

for j = 1:length(nuValueNames)
    
    outputValues = ['Eb_1_beta_0p5_nu_', nuValueNames{j}, '_kf_0p01_L0_0p125_homoggrowth'];
        
    load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
    % load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
    % load([outputDirectory, 'stresses_', outputValues,'.mat'], 'stressSols') % Foundation stresses
    load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
    load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
    
    index = find(times == timeSample, 1);
    
    % Initialise the current solution, gamma and the stresses
    initSol.x = Sols{index}(1,:);
    initSol.y = Sols{index}(2:(end - 1),:);
    initP = Sols{index}(end,:);
    
    % initGamma = gammaSols{end - 1};
    initGamma = parameters.gamma;
    % initP = stressSols{index};
    
    % parameters.gamma = initGamma(2, 1);
    parameters.gamma = 1 + timeSample*parameters.g;
    parameters.P = initP;
    
    solMesh = initSol.x;
    
    % Set the solver options for bvp4c
    maxPoints = 1e6;
    
    % Set the tolerances and max. number of mesh points
    solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');
    % Update the stresses until self-intersection or they tend to a steady state
    
    TMax = 2.0;
    dt = parameters.dt;
    newTimes = timeSample:dt:TMax;
    numSols = length(newTimes);
    
    SolsWithoutGrowth = cell(numSols, 1);
    
    solOld = initSol;
    
    SolsWithoutGrowth{1} = Sols{index};
    
    tic
    
    for i = 2:numSols
    
    % Update the solution
    [solMeshNew, solNew, PNew] = UpdateStandardLinearSolidSolutionWithoutGrowth(solMesh, solOld, parameters, solOptions);
    
    % Stop the solution the curve self-intersects or the stress tends to a
    % steady state
    if ( (HasRodHitSelfContact(solNew, parameters))||(norm(PNew - parameters.P, 2) < 1e-5) )
        
        SolsWithoutGrowth = SolsWithoutGrowth(1:(i - 1));
        newTimes = newTimes(1:(i - 1));
        
        break
    end
    
    % Update the solutions and the spring stresses
    parameters.P = interp1(solOld.x, PNew, solNew.x);
    solOld = solNew;
    
    SolsWithoutGrowth{i} = [solOld.x; solOld.y; parameters.P];
    
    solMesh = solMeshNew;
    
    end
    
    toc
    
    outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/StandardLinearSolid/';
    outputValues = ['Eb_1_beta_0p5_nu_', nuValueNames{j}, '_kf_0p01_L0_0p125_relaxation_homoggrowth'];
    save([outputDirectory, 'sols_', outputValues, '.mat'], 'SolsWithoutGrowth') % Solutions
%     save([outputDirectory, 'stresses_', outputValues,'.mat'], 'stressSolsWithoutGrowth') % Foundation stresses
    save([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
    
end

%%
outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/StandardLinearSolid/';
outputValues = ['Eb_1_beta_0p5_nu_1_kf_0p01_L0_0p125_relaxation_homoggrowth'];
load([outputDirectory, 'sols_', outputValues, '.mat'], 'SolsWithoutGrowth') % Solutions
%     save([outputDirectory, 'stresses_', outputValues,'.mat'], 'stressSolsWithoutGrowth') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times

for i = [1, length(newTimes)]
    
    hold on
    plot(SolsWithoutGrowth{i}(3,:), SolsWithoutGrowth{i}(4,:))
    
end
