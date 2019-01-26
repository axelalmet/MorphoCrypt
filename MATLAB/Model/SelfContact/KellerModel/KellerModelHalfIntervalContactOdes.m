function Odes = KellerModelHalfIntervalContactOdes(x, M, region, modelParams)

P = modelParams.P;
L = modelParams.L;

    function dMdS = ContactEqns(x, M, region)
        
        % Define the state variables
        S = M(1,:);
        x = M(2,:);
        y = M(3,:);
        Q = M(4,:);
        N = M(5,:);
        theta = M(6,:);
        k = M(7,:);
        sc = M(8,:);
        fc = M(9,:);
        
        dscdS = zeros(1, length(S));
        dfcdS = zeros(1, length(S));
        
        switch region
            case 1
                
                dSdS = (sc).*ones(1, length(S));
                dxdS = (sc).*cos(theta);
                dydS = (sc).*sin(theta);
                dQdS = -(sc).*k.*N;
                dNdS = (sc).*(k.*Q + P);
                dthetadS = (sc).*k;
                dkdS = (sc).*N;
                
            case 2
                
                dSdS = (L - sc).*ones(1, length(S));
                dxdS = (L - sc).*cos(theta);
                dydS = (L - sc).*sin(theta);
                dQdS = -(L - sc).*k.*N;
                dNdS = (L - sc).*(k.*Q + P);
                dthetadS = (L - sc).*k;
                dkdS = (L - sc).*N;    
                
            otherwise
                
                error('MATLAB:contactOdes:BadRegionIndex','Incorrect region index: %d',region);

        end
        
        dMdS = [dSdS; dxdS; dydS; dQdS; dNdS; dthetadS; dkdS; dscdS; dfcdS];
    end

Odes = ContactEqns(x, M, region);

end
