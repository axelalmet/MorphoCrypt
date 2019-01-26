function Odes = LinearElasticFoundationContactOdes(x, M, region, MOld, modelParams)

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
        A = M(8,:);
        sc = M(9,:);
        pC = M(10,:);
%         fC = M(11,:);

        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gammaInterp = interparc(x/2, MOld.x, gamma, 'linear');
            gamma = gammaInterp(:, 2)';
        end
        
        if (length(K) > 1)
            KInterp = interparc(x/2, MOld.x, K, 'linear');
            K = KInterp(:,2)';
        end
        
        if (length(Es) > 1)
            EsInterp = interparc(x/2, MOld.x, Es, 'linear');
            Es = EsInterp(:,2)';
        end
        
        if (length(Eb) > 1)
            EbInterp = interparc(x/2, MOld.x, Eb, 'linear');
            Eb = EbInterp(:,2)';
        end
        
        % If the model is extensible, set alpha to the tension, otherwise,
        % set alpha = 1.
        if (ext == 1)
            alpha = 1 + (F.*cos(theta) + G.*sin(theta))./Es;
        else
            alpha = 1;
        end
        
        dsCdS = zeros(1, length(S));
        dpCdS = zeros(1, length(S));
%         dfCdS = zeros(1, length(S));
        
        switch region
            case 1 % ODEs for [0, sc1]
 
                dSdS = (sc).*ones(1, length(S));
                dxdS = (sc).*gamma.*alpha.*cos(theta);
                dydS = (sc).*gamma.*alpha.*sin(theta);
                dFdS = (sc).*zeros(1, length(S));
                dGdS = (sc).*K.*Y;
                dthetadS = (sc).*gamma.*m./Eb;
                dmdS = (sc).*gamma.*alpha.*(F.*sin(theta) - G.*cos(theta));
                dAdS = zeros(1, length(S));

            case 2 % ODEs for [sc1, sc2]
                
                dSdS = (L - sc).*ones(1, length(S));
                dxdS = (L - sc).*gamma.*alpha.*cos(theta);
                dydS = (L - sc).*gamma.*alpha.*sin(theta);
                dFdS = -(L - sc).*pC.*dydS;
                dGdS = (L - sc).*(K.*Y + pC.*dxdS);
%                 dFdS = (L - sc).*zeros(1, length(S));
%                 dGdS = (L - sc).*K.*Y;
                dthetadS = (L - sc).*gamma.*m./Eb;
                dmdS = (L - sc).*gamma.*alpha.*(F.*sin(theta) - G.*cos(theta));
                dAdS = (L - sc).*X.*dydS;
                
            otherwise
                error('MATLAB:contactOdes:BadRegionIndex','Incorrect region index: %d',region);
        end
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS; dAdS; dsCdS; dpCdS];
%         dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS; dAdS; dsCdS; dpCdS; dfCdS];
    end

Odes = ContactEqns(x, M, region);

end
