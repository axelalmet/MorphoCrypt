function solMeshNew = LocallyRefineMesh(solMeshOld, currentSol)
% Function to locally refine mesh in a rather simplified manner, i.e. just
% refine the mesh where the solution changes drastically.

solMeshNew = solMeshOld; % Initialise the new mesh

for i = 2:length(solMeshOld)
    
    if (abs(currentSol(i) - currentSol(i - 1)) > 1.5)
        
        indexOne = find(solMeshNew == solMeshOld(i - 1), 1);
        indexTwo = find(solMeshNew == solMeshOld(i), 1);
        
        solMeshNew = [solMeshNew(1:(indexOne - 1)), ...
                    linspace(solMeshNew(indexOne), solMeshNew(indexTwo), 2), ...
                    solMeshNew((indexTwo + 1):end)];
    end
    
end
