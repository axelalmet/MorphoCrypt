function Odes = RemodellingFoundationSelfContactAlongRegionMiddleOdes(x, M, MOld, modelParams, sc1, sc2)

K = modelParams.K;
gamma = modelParams.gamma;

contIndex = find(MOld.x == 1);
SOld = MOld.y(1, [1:contIndex, (contIndex + 2):end]);

Py = modelParams.Py;

    function dMdS = InnerContactEqns(x, M)
        
        % Define the state variables
        S = M(1,:);
        Y = M(3,:);
        
        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gammaOld = gamma([1:contIndex, (contIndex + 2):end]);
            gamma = interp1(SOld, gammaOld, S);
        end
        
        if (length(Py) > 1)
            PyOld = Py([1:contIndex, (contIndex + 2):end]);
            Py = interp1(SOld, PyOld, S);
        end
        
        l = sc2 - sc1;
        
        dSdS = (l).*ones(1, length(S));
        dxdS = zeros(1, length(S));
        dydS = -(l).*gamma.*ones(1, length(S));
        dFdS = zeros(1, length(S));
        dGdS = (l).*gamma.*K.*(Y - Py);
        dthetadS = zeros(1, length(S));
        dmdS = zeros(1, length(S));
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = InnerContactEqns(x, M);

end
