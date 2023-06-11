function [t,states,algebraic] = INa(x)
    global num_alg_var;
    
    % Initialize constants and state variables
    states(:,1)
    states(:,2)
    STATES(:,20) = 0.713483e-6;  % ONa; Open state of fast Na+ channel

    constants(:,56) = 13;  % GNa; Maximun fast Na+ current conductance:mS/uF

end

function [rates,algebraic] = compute_rates(x,t,states)
    % Solve differential equations
    global num_alg_var;
    states_size = size(states);
    states_num_col = states_size(2);
    if ()

    else
        states_num_rows = states_size(1);
    
        % A51; alpha_Na11
        ALGEBRAIC(:,14) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./17.0000)+ 0.200000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));

        % A43; CNa2 (C2)
        RATES(:,22) = ( ALGEBRAIC(:,14).*ALGEBRAIC(:,4)+ ALGEBRAIC(:,38).*STATES(:,21)+ ALGEBRAIC(:,42).*STATES(:,26)) - ( ALGEBRAIC(:,36).*STATES(:,22)+ ALGEBRAIC(:,27).*STATES(:,22)+ ALGEBRAIC(:,44).*STATES(:,22));
        % A44; CNa1 (C1)
        RATES(:,21) = ( ALGEBRAIC(:,27).*STATES(:,22)+ ALGEBRAIC(:,40).*STATES(:,20)+ ALGEBRAIC(:,42).*STATES(:,25)) - ( ALGEBRAIC(:,38).*STATES(:,21)+ ALGEBRAIC(:,32).*STATES(:,21)+ ALGEBRAIC(:,44).*STATES(:,21));
        % A45; ONa (O)
        RATES(:,20) = ( ALGEBRAIC(:,32).*STATES(:,21)+ ALGEBRAIC(:,48).*STATES(:,25)) - ( ALGEBRAIC(:,40).*STATES(:,20)+ ALGEBRAIC(:,46).*STATES(:,20));
        % A46; IFNa (IF)
        RATES(:,25) = ( ALGEBRAIC(:,46).*STATES(:,20)+ ALGEBRAIC(:,44).*STATES(:,21)+ ALGEBRAIC(:,52).*STATES(:,23)+ ALGEBRAIC(:,27).*STATES(:,26)) - ( ALGEBRAIC(:,48).*STATES(:,25)+ ALGEBRAIC(:,42).*STATES(:,25)+ ALGEBRAIC(:,50).*STATES(:,25)+ ALGEBRAIC(:,38).*STATES(:,25));
        % A47; I1Na (I1)
        RATES(:,23) = ( ALGEBRAIC(:,50).*STATES(:,25)+ ALGEBRAIC(:,56).*STATES(:,24)) - ( ALGEBRAIC(:,52).*STATES(:,23)+ ALGEBRAIC(:,54).*STATES(:,23));
        % A48; I2Na (I2)
        RATES(:,24) =  ALGEBRAIC(:,54).*STATES(:,23) -  ALGEBRAIC(:,56).*STATES(:,24);
        % A49; ICNa2 (IC2)
        RATES(:,26) = ( ALGEBRAIC(:,14).*STATES(:,27)+ ALGEBRAIC(:,38).*STATES(:,25)+ ALGEBRAIC(:,44).*STATES(:,22)) - ( ALGEBRAIC(:,36).*STATES(:,26)+ ALGEBRAIC(:,27).*STATES(:,26)+ ALGEBRAIC(:,42).*STATES(:,26));
        % A50; ICNa3 (IC3)
        RATES(:,27) = ( ALGEBRAIC(:,36).*STATES(:,26)+ ALGEBRAIC(:,44).*ALGEBRAIC(:,4)) - ( ALGEBRAIC(:,14).*STATES(:,27)+ ALGEBRAIC(:,42).*STATES(:,27));

        % A42; CNa3 (C3)
        ALGEBRAIC(:,4) = 1.00000 - (STATES(:,20)+STATES(:,21)+STATES(:,22)+STATES(:,25)+STATES(:,23)+STATES(:,24)+STATES(:,26)+STATES(:,27));

        % A52; alpha_Na12
        ALGEBRAIC(:,27) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./15.0000)+ 0.230000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
        % A53; alpha_Na13
        ALGEBRAIC(:,32) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./12.0000)+ 0.250000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
        % A54; beta_Na11
        ALGEBRAIC(:,36) =  0.191700.*exp( - (ALGEBRAIC(:,72)+2.50000)./20.3000);
        % A55; beta_Na12
        ALGEBRAIC(:,38) =  0.200000.*exp( - (ALGEBRAIC(:,72) - 2.50000)./20.3000);
        % A56; beta_Na13
        ALGEBRAIC(:,40) =  0.220000.*exp( - (ALGEBRAIC(:,72) - 7.50000)./20.3000);
        % A57; alpha_Na3
        ALGEBRAIC(:,42) =  7.00000e-07.*exp( - (ALGEBRAIC(:,72)+7.00000)./7.70000);
        % A58; beta_Na3
        ALGEBRAIC(:,44) = 0.00840000+ 2.00000e-05.*(ALGEBRAIC(:,72)+7.00000);
        % A59; alpha_Na2
        ALGEBRAIC(:,46) = 1.00000./( 0.188495.*exp( - (ALGEBRAIC(:,72)+7.00000)./16.6000)+0.393956);
        % A60; beta_Na2
        ALGEBRAIC(:,48) = ( ALGEBRAIC(:,32).*ALGEBRAIC(:,46).*ALGEBRAIC(:,42))./( ALGEBRAIC(:,40).*ALGEBRAIC(:,44));
        % A61; alpha_Na4
        ALGEBRAIC(:,50) = ALGEBRAIC(:,46)./1000.00;
        % A62; beta_Na4
        ALGEBRAIC(:,52) = ALGEBRAIC(:,42);
        % A63; alpha_Na5
        ALGEBRAIC(:,54) = ALGEBRAIC(:,46)./95000.0;
        % A64; beta_Na5
        ALGEBRAIC(:,56) = ALGEBRAIC(:,42)./50.0000;        
    end
end

function algebraic = compute_alg(x,algebraic,states)
    % Compute algebraic equations related to INa

end

function vc = volt_clamp()
    % Generate voltage-clamp protocol
end

