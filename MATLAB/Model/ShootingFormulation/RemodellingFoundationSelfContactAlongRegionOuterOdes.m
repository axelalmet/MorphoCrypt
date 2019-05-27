function Odes = RemodellingFoundationSelfContactAlongRegionOuterOdes(x, M, MOld, modelParams, sc2)

K = modelParams.K;
gamma = modelParams.gamma;
Eb = modelParams.Eb;

contIndexOne = find(MOld.x == 1, 1);
contIndexTwo = find(MOld.x == 2, 1);
SOld = MOld.y(1, [1:contIndex, (contIndex + 2):contIndexTwo, (contIndexTwo + 2):end]);

Px = modelParams.Px;
Py = modelParams.Py;

    function dMdS = InnerContactEqns(x, M)
        
        % Define the state variables
        S = M(1,:);
        X = M(2,:);
        Y = M(3,:);
        F = M(4,:);
        G = M(5,:);
        theta = M(6,:);
        m = M(7,:);
        
        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gammaOld = gamma([1:contIndex, (contIndex + 2):end]);
            gamma = interp1(SOld, gammaOld, S);
        end
        
        if (length(Eb) > 1)
            EbOld = Eb([1:contIndex, (contIndex + 2):end]);
            Eb = interp1(SOld, EbOld, S);
        end
        
        if (length(Px) > 1)
            PxOld = Px([1:contIndex, (contIndex + 2):end]);
            Px = interp1(SOld, PxOld, S);
        end
        
        if (length(Py) > 1)
            PyOld = Py([1:contIndex, (contIndex + 2):end]);
            Py = interp1(SOld, PyOld, S);
        end
        
        l = L - sc2;
        
        dSdS = (l).*ones(1, length(S));
        dxdS = (l).*gamma.*cos(theta);
        dydS = (l).*gamma.*sin(theta);
        dFdS = (l).*gamma.*K.*(X - Px);
        dGdS = (l).*gamma.*K.*(Y - Py);
        dthetadS = (l).*gamma.*(m./Eb);
        dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = InnerContactEqns(x, M);

end
