function Odes = HalfIntervalOdes(x, M, region, MOld, modelParams)

K = modelParams.K;
L = modelParams.L;
gamma = modelParams.gamma;
Es = modelParams.Es;
Eb = modelParams.Eb;
ext = modelParams.ext;

    function dMdS = ContactEqns(x, M, region)
        
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
            gammaInterp = interparc(x, MOld.x, gamma, 'linear');
            gamma = gammaInterp(:, 2)';
        end
        
        if (length(K) > 1)
            KInterp = interparc(x, MOld.x, K, 'linear');
            K = KInterp(:,2)';
        end
        
        if (length(Es) > 1)
            EsInterp = interparc(x, MOld.x, Es, 'linear');
            Es = EsInterp(:,2)';
        end
        
        if (length(Eb) > 1)
            EbInterp = interparc(x, MOld.x, Eb, 'linear');
            Eb = EbInterp(:,2)';
        end
        
        % If the model is extensible, set alpha to the tension, otherwise,
        % set alpha = 1.
        if (ext == 1)
            alpha = 1 + (F.*cos(theta) + G.*sin(theta))./Es;
        else
            alpha = 1;
        end
        
        dSdS = L.*ones(1, length(S));
        dxdS = L.*gamma.*alpha.*cos(theta);
        dydS = L.*gamma.*alpha.*sin(theta);
        dFdS = L.*K.*(X - S);
        dGdS = L.*K.*Y;
        dthetadS = L.*gamma.*m./Eb;
        dmdS = L.*gamma.*alpha.*(F.*sin(theta) - G.*cos(theta));
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
    end

Odes = ContactEqns(x, M, region);

end
