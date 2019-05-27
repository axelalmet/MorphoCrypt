function Odes = UnbuckledRodRemodellingFoundationSloughingOdes(x, M, region, MOld, parameters)

K = parameters.K;
L = parameters.L;
gamma = parameters.gamma;

Px = parameters.Px;

contIndex = find(MOld.x == 1, 1);
SOld = MOld.y(1, [1:contIndex, (contIndex + 2):end]);


    function dMdS = FoundationEqns(x, M, region)
        
        % Define the state variables
        S = M(1,:);
        X = M(2,:);
        Y = M(3,:);
        F = M(4,:);
        
        if (length(gamma) > 1)
            gammaOld = gamma([1:contIndex, (contIndex + 2):end]);
            gamma = interp1(SOld, gammaOld, S);
        end
        
        if (length(K) > 1)
            KOld = K([1:contIndex, (contIndex + 2):end]);
            K = interp1(SOld, KOld, S);
        end
        
        if (length(Px) > 1)
            PxOld = Px([1:contIndex, (contIndex + 2):end]);
            Px = interp1(SOld, PxOld, S);
        end
        
        switch region
            
            case 1
                
                l = L;
                
                dSdS = l.*ones(1, length(S));
                dxdS = l.*gamma.*ones(1, length(S));
                dydS = zeros(1, length(S));
                dFdS = l.*K.*gamma.*(X - Px);
                
            case 2
                
                l = 0.5 - L;
                
                dSdS = l.*ones(1, length(S));
                dxdS = l.*gamma.*ones(1, length(S));
                dydS = zeros(1, length(S));
                dFdS = zeros(1, length(S));
                
            otherwise
                
                error('MATLAB:sloughingodes:BadRegionIndex','Incorrect region index: %d',region);
                
        end
        
        dMdS = [dSdS; dxdS; dydS; dFdS];
        
    end

Odes = FoundationEqns(x, M, region);

end
