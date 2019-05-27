function Odes = RemodellingFoundationSelfContactAtPointOdes(x, M, MOld, modelParams, sc, fc)

K = modelParams.K;
gamma = modelParams.gamma;
Eb = modelParams.Eb;

contIndex = find(MOld.x == 1, 1);
xOld = MOld.x(1, [1:contIndex, (contIndex + 2):end]);

Px = modelParams.Px;
Py = modelParams.Py;

    function dMdS = InnerContactEqns(x, M)
        
        % Define the state variables
        X = M(1,:);
        Y = M(2,:);
        F = M(3,:);
        G = M(4,:);
        theta = M(5,:);
        m = M(6,:);
        
        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gammaOld = gamma([1:contIndex, (contIndex + 2):end]);
            gamma = interp1(xOld, gammaOld, x);
        end
        
        if (length(Eb) > 1)
            EbOld = Eb([1:contIndex, (contIndex + 2):end]);
            Eb = interp1(xOld, EbOld, x);
        end
        
        if (length(Px) > 1)
            PxOld = Px([1:contIndex, (contIndex + 2):end]);
            Px = interp1(xOld, PxOld, x);
        end
        
        if (length(Py) > 1)
            PyOld = Py([1:contIndex, (contIndex + 2):end]);
            Py = interp1(xOld, PyOld, x);
        end
        
                
        dxdS = gamma.*cos(theta);
        dydS = gamma.*sin(theta);
        dFdS = gamma.*K.*(X - Px) - fc.*dirac(x - sc);
        dGdS = gamma.*K.*(Y - Py);
        dthetadS = gamma.*(m./Eb);
        dmdS = gamma.*(F.*sin(theta) - G.*cos(theta));
        
        dMdS = [dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = InnerContactEqns(x, M);

end
