function curveNew = InterpolateToNewMesh(meshOld, curveOld, meshNew)
% This function is designed to interpolate the curves to the meshes that
% are needed for bvp4c's multipoint functionality. For piecewise-continuous
% meshes, interp1 doesn't quite work, so we will split it up

% Find where the old and new meshes split
splitIndexOld = find(meshOld == 1, 1);
splitIndexNew = find(meshNew == 1, 1);

curveInterpOne = interp1(meshOld(1:splitIndexOld), curveOld(1:splitIndexOld), ...
                    meshNew(1:splitIndexNew));
                              
curveInterpTwo = interp1(meshOld((splitIndexOld + 1):end), curveOld((splitIndexOld + 1):end), ...
                    meshNew((splitIndexNew + 1):end)); 
         
curveNew = [curveInterpOne, curveInterpTwo];



