function Odes = RemodellingFoundationWithSloughingOdes(x, M, region, MOld, parameters)

K = parameters.K;
L = parameters.L;
gamma = parameters.gamma;
Eb = parameters.Eb;

Px = parameters.Px;
Py = parameters.Py;

contIndex = find(MOld.x == 1, 1);
SOld = MOld.y(1, [1:contIndex, (contIndex + 2):end]);


    function dMdS = FoundationEqns(x, M, region)
        
        % Define the state variables
        S = M(1,:);
        X = M(2,:);
        Y = M(3,:);
        F = M(4,:);
        G = M(5,:);
        theta = M(6,:);
        m = M(7,:);
        
        if (length(gamma) > 1)
            gammaOld = gamma([1:contIndex, (contIndex + 2):end]);
            gamma = interp1(SOld, gammaOld, S);
        end
        
        if (length(K) > 1)
            KOld = K([1:contIndex, (contIndex + 2):end]);
            K = interp1(SOld, KOld, S);
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
        
        switch region
            
            case 1
                
                l = L;
                
                dSdS = l.*ones(1, length(S));
                dxdS = l.*gamma.*cos(theta);
                dydS = l.*gamma.*sin(theta);
                dFdS = l.*K.*gamma.*(X - Px);
                dGdS = l.*K.*gamma.*(Y - Py);
                dthetadS = l.*gamma.*(m./Eb);
                dmdS = l.*gamma.*(F.*sin(theta) - G.*cos(theta));
                
            case 2
                
                l = 0.5 - L;
                
                dSdS = l.*ones(1, length(S));
                dxdS = l.*gamma.*ones(1, length(S));
                dydS = zeros(1, length(S));
                dFdS = zeros(1, length(S));
                dGdS = zeros(1, length(S));
                dthetadS = zeros(1, length(S));
                dmdS = zeros(1, length(S));
            otherwise
                
                error('MATLAB:sloughingodes:BadRegionIndex','Incorrect region index: %d',region);
                
        end
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = FoundationEqns(x, M, region);

end
