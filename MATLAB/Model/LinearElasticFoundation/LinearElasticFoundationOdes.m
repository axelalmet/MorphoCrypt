function Odes = LinearElasticFoundationOdes(x, M, MOld, parameters)

K = parameters.K;
L = parameters.L;
gamma = parameters.gamma;
Eb = parameters.Eb;
y0 = parameters.y0;
uHat = parameters.uHat;

    function dMdS = TensionBasedGrowthEqns(x, M)
        
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
            gamma = interp1(MOld.x, gamma, x);
        end
        
        if (length(K) > 1)
            K = interp1(MOld.x, K, x);
        end
        
        if (length(Eb) > 1)
            Eb = interp1(MOld.x, Eb, x);
        end
            
        if (length(uHat) > 1)
            uHat = interp1(MOld.x, uHat, x);
        end
        
        dSdS = L.*ones(1, length(S));
        dxdS = L.*gamma.*cos(theta);
        dydS = L.*gamma.*sin(theta);
        dFdS = L.*K.*(X - S);
        dGdS = L.*K.*Y;
        dthetadS = L.*gamma.*(m./Eb + uHat);
        dmdS = L.*gamma.*(F.*sin(theta) - G.*cos(theta));

        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
    end

Odes = TensionBasedGrowthEqns(x, M);

end
