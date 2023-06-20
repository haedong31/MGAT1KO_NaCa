clc
clearvars

[t,s,a] = INa();

figure('Color','w')
plot(t,a(:,16))

figure('Color','w')
plot(t,s)
legend(["ONa","CNa1","CNa2","I1Na","I2Na","IFNa","ICNa2","ICNa3"])

figure('Color','w')
plot(t,a)
legend(["alpha_Na11","alpha_Na12","alpha_Na13","beta_Na11","beta_Na12", ...
    "beta_Na13","alpha_Na3","beta_Na3","alpha_Na2","beta_Na2","alpha_Na4", ...
    "beta_Na4","alpha_Na5","beta_Na5","CNa3 (C3)","I_Na","Volt"])

function [t,states,algebraic] = INa()
    % Initialize state variables
    init_states = [];
    init_states(:,1) = 0.713483e-6; % ONa; Open state of fast Na+ channel
    init_states(:,2) = 0.279132e-3; % CNa1; Closed state of fast Na+ channel
    init_states(:,3) = 0.020752; % CNa2; Closed state of fast Na+ channel
    init_states(:,4) = 0.673345e-6; % I1Na; Slow inactivated state 1 of fast Na+ channel
    init_states(:,5) = 0.155787e-8; % I2Na; Slow inactivated state 2 of fast Na+ channel
    init_states(:,6) = 0.153176e-3; % IFNa; Fast inactivated state of fast Na+ channel
    init_states(:,7) = 0.0113879; % ICNa2; Cloesd-inactivated state of fast Na+ channel
    init_states(:,8) = 0.34278; % ICNa3; Cloesd-inactivated state of fast Na+ channel

    % Constant variables
    constants = [];
    constants(1) = 13.0; % GNa; Maximun fast Na+ current conductance :mS/uF
    constants(2) = 39.4843; % A41; ENa; 39.4843
    constants(3) = 5; % Nernst potential

    % Set options for ODE solver
    tspan = [0,150];
    options = odeset('RelTol',1e-06,'AbsTol',1e-06,'MaxStep',1);

    % Solve model using ODE solver
    [t,states] = ode15s(@(t,states)compute_rates(t,states,constants),tspan,init_states,options);

    % Compute algebraic variables
    [~,algebraic] = compute_rates(t,states,constants);
    algebraic = compute_alg(t,algebraic,states,constants);
end

function [rates,algebraic] = compute_rates(t,states,constants)
    % Solve differential equations
    % number of state: 8
    % number of algebraic variables: 17
    num_alg_var = 17;
    states_size = size(states);
    states_num_col = states_size(2);
    if (states_num_col == 1)
        states = states';
        algebraic = zeros(1,num_alg_var);
    else
        states_num_rows = states_size(1);
        algebraic = zeros(states_num_rows,num_alg_var);
        rates = zeros(states_num_rows,states_num_col);
    end

    % Stimulation voltage
    algebraic(:,17) = arrayfun(@(t)volt_clamp(t,constants(3)),t);
    
    % A51; alpha_Na11
    algebraic(:,1) = 3.80200./( 0.1027.*exp(-(algebraic(:,17)+2.50000)./17.0)+ 0.2.*exp(-(algebraic(:,17)+2.50000)./150.000));
    % A52; alpha_Na12
    algebraic(:,2) = 3.80200./( 0.1027.*exp(-(algebraic(:,17)+2.50000)./15.0)+ 0.23.*exp(-(algebraic(:,17)+2.50000)./150.000));
    % A53; alpha_Na13
    algebraic(:,3) = 3.80200./( 0.1027.*exp(-(algebraic(:,17)+2.50000)./12.0)+ 0.25.*exp(-(algebraic(:,17)+2.50000)./150.000));
    % A54; beta_Na11
    algebraic(:,4) =  0.1917.*exp(-(algebraic(:,17)+2.50000)./20.3000);
    % A55; beta_Na12
    algebraic(:,5) =  0.2.*exp(-(algebraic(:,17)-2.50000)./20.3000);
    % A56; beta_Na13
    algebraic(:,6) =  0.22.*exp(-(algebraic(:,17)-7.50000)./20.3000);
    % A57; alpha_Na3
    algebraic(:,7) =  7.0e-07.*exp(-(algebraic(:,17)+7.0)./7.7);
    % A58; beta_Na3
    algebraic(:,8) = 0.0084+ 2.0e-05.*(algebraic(:,17)+7.0);
    % A59; alpha_Na2
    algebraic(:,9) = 1.0./( 0.188495.*exp(-(algebraic(:,17)+7.0)./16.6000)+0.393956);
    % A60; beta_Na2
    algebraic(:,10) = ( algebraic(:,3).*algebraic(:,9).*algebraic(:,7))./( algebraic(:,6).*algebraic(:,8));
    % A61; alpha_Na4
    algebraic(:,11) = algebraic(:,9)./1000.00;
    % A62; beta_Na4
    algebraic(:,12) = algebraic(:,7);
    % A63; alpha_Na5
    algebraic(:,13) = algebraic(:,9)./95000.0;
    % A64; beta_Na5
    algebraic(:,14) = algebraic(:,7)./50.0000;

    % A42; CNa3 (C3); 176
    algebraic(:,15) = 1.0 - (states(:,1)+states(:,2)+states(:,3)+states(:,6)+states(:,4)+states(:,5)+states(:,7)+states(:,8));
    % A43; RATES(:,22); CNa2 (C2); 190
    rates(:,3) = (algebraic(:,1).*algebraic(:,15)+ algebraic(:,5).*states(:,2)+ algebraic(:,7).*states(:,7)) - (algebraic(:,4).*states(:,3)+ algebraic(:,2).*states(:,3)+ algebraic(:,8).*states(:,3));
    % A44; RATES(:,21); CNa1 (C1); 196
    rates(:,2) = (algebraic(:,2).*states(:,3)+ algebraic(:,6).*states(:,1)+ algebraic(:,7).*states(:,6)) - (algebraic(:,5).*states(:,2)+ algebraic(:,3).*states(:,2)+ algebraic(:,8).*states(:,2));
    % A49; RATES(:,26); ICNa2 (IC2); 198
    rates(:,7) = (algebraic(:,1).*states(:,8)+ algebraic(:,5).*states(:,6)+ algebraic(:,8).*states(:,3)) - (algebraic(:,4).*states(:,7)+ algebraic(:,2).*states(:,7)+ algebraic(:,7).*states(:,7));
    % A50; RATES(:,27); ICNa3 (IC3); 200
    rates(:,8) = (algebraic(:,4).*states(:,7)+ algebraic(:,8).*algebraic(:,15)) - (algebraic(:,1).*states(:,8)+ algebraic(:,7).*states(:,8));    
    % A45; RATES(:,20); ONa (O); 216
    rates(:,1) = (algebraic(:,3).*states(:,2)+ algebraic(:,10).*states(:,6)) - (algebraic(:,6).*states(:,1)+ algebraic(:,9).*states(:,1));
    % A46; RATES(:,25); nIFNa (IF); 222
    rates(:,6) = (algebraic(:,9).*states(:,1)+ algebraic(:,8).*states(:,2)+ algebraic(:,12).*states(:,4)+ algebraic(:,2).*states(:,7)) - (algebraic(:,10).*states(:,6)+ algebraic(:,7).*states(:,6)+ algebraic(:,11).*states(:,6)+ algebraic(:,5).*states(:,6));
    % A47; RATES(:,23); I1Na (I1); 242
    rates(:,4) = (algebraic(:,11).*states(:,6)+ algebraic(:,14).*states(:,5)) - (algebraic(:,12).*states(:,4)+ algebraic(:,13).*states(:,4));
    % A48; RATES(:,24); I2Na (I2); 244
    rates(:,5) =  algebraic(:,13).*states(:,4) -  algebraic(:,14).*states(:,5);
    % A40; I_Na
    algebraic(:,16) =  constants(1).*states(:,1).*(algebraic(:,17) - constants(2));   
    
    rates = rates';
end

function algebraic = compute_alg(t,algebraic,states,constants)
    % Compute algebraic equations related to INa
    algebraic(:,17) = arrayfun(@(t)volt_clamp(t,constants(3)),t);
    
    % 1. alpha11; A51; ALGEBRAIC(:,14)
    algebraic(:,1) = 3.802./( 0.1027.*exp(-(algebraic(:,17)+x(1))./17.0)+ 0.2.*exp(-(algebraic(:,17)+x(1))./150.0));
    % 2. beta11; A54; ALGEBRAIC(:,36)
    algebraic(:,4) =  0.1917.*exp(-(algebraic(:,17)+x(2))./20.3);
    % 3. alpha12; A52; ALGEBRAIC(:,27)
    algebraic(:,2) = 3.802./( 0.1027.*exp(-(algebraic(:,17)+x(3))./15.0)+ 0.23.*exp(-(algebraic(:,17)+x(3))./150.0));
    % 4. beta12; A55; ALGEBRAIC(:,38)
    algebraic(:,5) =  0.2.*exp(-(algebraic(:,17)-x(4))./20.3);
    % 5. alpha13; A53; ALGEBRAIC(:,32)
    algebraic(:,3) = 3.802./( 0.1027.*exp(-(algebraic(:,17)+x(5))./12.0)+ 0.25.*exp(-(algebraic(:,17)+x(5))./150.0));
    % 6. beta13; A56; ALGEBRAIC(:,40)
    algebraic(:,6) =  0.22.*exp(-(algebraic(:,17)-x(6))./20.3);
    % 7. alpha111
    3.802./(0.1027.*exp(-(algebraic(:,17)+x(7))./17.0) + 0.2.*exp(-(algebraic(:,17)+x(7))./150.0));
    % 8. alpha112
    3.802./(0.1027.*exp(-(algebraic(:,17)+x(8))./15.0) + 0.23.*exp(-(algebraic(:,17)+x(8))./150.0));
    % 9. beta111
    0.1917.*exp(-(algebraic(:,17)+x(9))./20.3);
    % 10. beta112
    0.2.*exp(-(algebraic(:,17)-x(10))./20.3);
    % 11. alpha31; A57; ALGEBRAIC(:,42)
    algebraic(:,7) = 7.0e-07.*exp(-(algebraic(:,17)+x(11))./x(23));
    % 12. beta31; A58; ALGEBRAIC(:,44)
    algebraic(:,8) = x(12)+ 2.0e-05.*(algebraic(:,17)+7.0);
    % 13. alpha32
    7.0e-07.*exp(-(algebraic(:,17)+x(13))./x(24));
    % 14. beta32
    x(14) + 2.0e-05.*(algebraic(:,17)+7.0);
    % 15. alpha33
    7.0e-07.*exp(-(algebraic(:,17)+x(15))./x(25));
    % 16. beta33
    x(16) + 2.0e-05.*(algebraic(:,17)+7.0);
    % 17. alpha2; A59; ALGEBRAIC(:,46)
    algebraic(:,9) = 1.0./( 0.188495.*exp(-(algebraic(:,17)+x(17))./x(22))+0.393956);
    % 18. beta2; A60; ALGEBRAIC(:,48)
    % alpha3; algebraic(:,7) -> alpha33 algebraic(:,10) = ( algebraic(:,3).*algebraic(:,9).*algebraic(:,7))./( algebraic(:,6).*algebraic(:,8));
    % 19. alpha4; A61; ALGEBRAIC(:,50)
    algebraic(:,11) = x(18).*algebraic(:,9);
    % 20. beta4; A62; ALGEBRAIC(:,52)
    algebraic(:,12) = x(19).*algebraic(:,7);
    % 21. alpha5; A63; ALGEBRAIC(:,54)
    algebraic(:,13) = x(20).*algebraic(:,9);
    % 22. beta5; A64; ALGEBRAIC(:,56)
    algebraic(:,14) = x(21).*algebraic(:,7);
    
    % A42; CNa3 (C3); ALGEBRAIC(:,4)
    algebraic(:,15) = 1.0 - (states(:,1)+states(:,2)+states(:,3)+states(:,6)+states(:,4)+states(:,5)+states(:,7)+states(:,8));
    % A40; I_Na; ALGEBRAIC(:,58)
    algebraic(:,16) =  constants(1).*states(:,1).*(algebraic(:,17) - constants(2));

    % x(18) = 1/1000.0; x(19) = 1.0; x(20) 1/95000.0; x(21) = 1/50.0;
end

function vc = volt_clamp(t,p)
    % Generate voltage-clamp protocol
    holdt = 5;
    pt = 120;
    holdp = -100;
    if (t < holdt)
        vc = holdp;
    elseif (t >= holdt) && (t <= pt)
        vc = p;
    else
        vc = holdp;
    end
end
