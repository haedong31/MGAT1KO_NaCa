clc
clearvars

x0 = [2.5,2.5,2.5,-2.5,2.5,-7.5,2.5,2.5,2.5,-2.5,7,0.0084,7,0.0084,7,0.0084,7,1/1000.0,1.0,1/95000.0,1/50.0,16.6,7.7,7.7,7.7];
[t,s,a] = INa(x0);

figure('Color','w')
plot(t,a(:,24))
title("I_{Na}")

figure('Color','w')
plot(t,s)
legend(["ONa","CNa1","CNa2","I1Na","I2Na","IFNa","ICNa2","ICNa3"])

% figure('Color','w')
% plot(t,a)
% legend(["alpha_Na11","alpha_Na12","alpha_Na13","beta_Na11","beta_Na12", ...
%     "beta_Na13","alpha_Na3","beta_Na3","alpha_Na2","beta_Na2","alpha_Na4", ...
%     "beta_Na4","alpha_Na5","beta_Na5","CNa3 (C3)","I_Na","Volt"])

function [t,states,algebraic] = INa(x)
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
    [t,states] = ode15s(@(t,states)compute_rates(t,states,constants,x),tspan,init_states,options);

    % Compute algebraic variables
    [~,algebraic] = compute_rates(t,states,constants,x);
    algebraic = compute_alg(t,algebraic,states,constants,x);
end

function [rates,alg] = compute_rates(t,states,constants,x)
    % Solve differential equations
    % number of state: 8
    % number of algebraic variables: 25
    num_alg_var = 25;
    states_size = size(states);
    states_num_col = states_size(2);
    if (states_num_col == 1)
        states = states';
        alg = zeros(1,num_alg_var);
    else
        states_num_rows = states_size(1);
        alg = zeros(states_num_rows,num_alg_var);
        rates = zeros(states_num_rows,states_num_col);
    end

    % Stimulation voltage
    alg(:,25) = arrayfun(@(t)volt_clamp(t,constants(3)),t);
    
    % 1. alpha11; A51; ALGEBRAIC(:,14)
    alg(:,1) = 3.802./( 0.1027.*exp(-(alg(:,25)+x(1))./17.0)+ 0.2.*exp(-(alg(:,25)+x(1))./150.0));
    % 2. beta11; A54; ALGEBRAIC(:,36)
    alg(:,2) =  0.1917.*exp(-(alg(:,25)+x(2))./20.3);
    % 3. alpha12; A52; ALGEBRAIC(:,27)
    alg(:,3) = 3.802./( 0.1027.*exp(-(alg(:,25)+x(3))./15.0)+ 0.23.*exp(-(alg(:,25)+x(3))./150.0));
    % 4. beta12; A55; ALGEBRAIC(:,38)
    alg(:,4) =  0.2.*exp(-(alg(:,25)-x(4))./20.3);
    % 5. alpha13; A53; ALGEBRAIC(:,32)
    alg(:,5) = 3.802./( 0.1027.*exp(-(alg(:,25)+x(5))./12.0)+ 0.25.*exp(-(alg(:,25)+x(5))./150.0));
    % 6. beta13; A56; ALGEBRAIC(:,40)
    alg(:,6) =  0.22.*exp(-(alg(:,25)-x(6))./20.3);
    % 7. alpha111
    alg(:,7) = 3.802./(0.1027.*exp(-(alg(:,25)+x(7))./17.0) + 0.2.*exp(-(alg(:,25)+x(7))./150.0));
    % 8. alpha112
    alg(:,8) = 3.802./(0.1027.*exp(-(alg(:,25)+x(8))./15.0) + 0.23.*exp(-(alg(:,25)+x(8))./150.0));
    % 9. beta111
    alg(:,9) = 0.1917.*exp(-(alg(:,25)+x(9))./20.3);
    % 10. beta112
    alg(:,10) = 0.2.*exp(-(alg(:,25)-x(10))./20.3);
    % 11. alpha31; A57; ALGEBRAIC(:,42)
    alg(:,11) = 7.0e-07.*exp(-(alg(:,25)+x(11))./x(23));
    % 12. beta31; A58; ALGEBRAIC(:,44)
    alg(:,12) = x(12)+ 2.0e-05.*(alg(:,25)+7.0);
    % 13. alpha32
    alg(:,13) = 7.0e-07.*exp(-(alg(:,25)+x(13))./x(24));
    % 14. beta32
    alg(:,14) = x(14) + 2.0e-05.*(alg(:,25)+7.0);
    % 15. alpha33
    alg(:,15) = 7.0e-07.*exp(-(alg(:,25)+x(15))./x(25));
    % 16. beta33
    alg(:,16) = x(16) + 2.0e-05.*(alg(:,25)+7.0);
    % 17. alpha2; A59; ALGEBRAIC(:,46)
    alg(:,17) = 1.0./( 0.188495.*exp(-(alg(:,25)+x(17))./x(22))+0.393956);
    % 18. beta2; A60; ALGEBRAIC(:,48) alpha13.*alpha2.*alpha33./(beta13.*beta33)
    alg(:,18) = (alg(:,5).*alg(:,17).*alg(:,15))./(alg(:,6).*alg(:,16));
    % 19. alpha4; A61; ALGEBRAIC(:,50)
    alg(:,19) = x(18).*alg(:,9);
    % 20. beta4; A62; ALGEBRAIC(:,52)
    alg(:,20) = x(19).*alg(:,11);
    % 21. alpha5; A63; ALGEBRAIC(:,54)
    alg(:,21) = x(20).*alg(:,9);
    % 22. beta5; A64; ALGEBRAIC(:,56)
    alg(:,22) = x(21).*alg(:,11);

    % A42; CNa3 (C3); 176
    alg(:,23) = 1.0 - (states(:,1)+states(:,2)+states(:,3)+states(:,6)+states(:,4)+states(:,5)+states(:,7)+states(:,8));
    % A43; RATES(:,22); CNa2 (C2); 190
    rates(:,3) = (alg(:,1).*alg(:,23)+alg(:,4).*states(:,2)+alg(:,13).*states(:,7)) - (alg(:,2).*states(:,3)+alg(:,3).*states(:,3)+alg(:,14).*states(:,3));
    % A44; RATES(:,21); CNa1 (C1); 196
    rates(:,2) = (alg(:,3).*states(:,3)+alg(:,6).*states(:,1)+alg(:,15).*states(:,6)) - (alg(:,4).*states(:,2)+alg(:,5).*states(:,2)+alg(:,16).*states(:,2));
    % A49; RATES(:,26); ICNa2 (IC2); 198
    rates(:,7) = (alg(:,7).*states(:,8)+alg(:,10).*states(:,6)+alg(:,14).*states(:,3)) - (alg(:,9).*states(:,7)+alg(:,8).*states(:,7)+ alg(:,13).*states(:,7));
    % A50; RATES(:,27); ICNa3 (IC3); 200
    rates(:,8) = (alg(:,9).*states(:,7)+alg(:,12).*alg(:,23)) - (alg(:,7).*states(:,8)+alg(:,11).*states(:,8));    
    % A45; RATES(:,20); ONa (O); 216
    rates(:,1) = (alg(:,5).*states(:,2)+alg(:,18).*states(:,6)) - (alg(:,6).*states(:,1)+alg(:,17).*states(:,1));
    % A46; RATES(:,25); nIFNa (IF); 222
    rates(:,6) = (alg(:,17).*states(:,1)+alg(:,16).*states(:,2)+alg(:,20).*states(:,4)+alg(:,8).*states(:,7)) - (alg(:,10).*states(:,6)+alg(:,15).*states(:,6)+alg(:,19).*states(:,6)+alg(:,18).*states(:,6));
    % A47; RATES(:,23); I1Na (I1); 242
    rates(:,4) = (alg(:,19).*states(:,6)+alg(:,22).*states(:,5)) - (alg(:,20).*states(:,4)+alg(:,21).*states(:,4));
    % A48; RATES(:,24); I2Na (I2); 244
    rates(:,5) =  alg(:,21).*states(:,4) - alg(:,22).*states(:,5);
    % A40; I_Na
    alg(:,24) = constants(1).*states(:,1).*(alg(:,25) - constants(2));   
    
    rates = rates';
end

function alg = compute_alg(t,alg,states,constants,x)
    % Compute algebraic equations related to INa
    alg(:,25) = arrayfun(@(t)volt_clamp(t,constants(3)),t);
    
    % 1. alpha11; A51; ALGEBRAIC(:,14)
    alg(:,1) = 3.802./( 0.1027.*exp(-(alg(:,25)+x(1))./17.0)+ 0.2.*exp(-(alg(:,25)+x(1))./150.0));
    % 2. beta11; A54; ALGEBRAIC(:,36)
    alg(:,2) =  0.1917.*exp(-(alg(:,25)+x(2))./20.3);
    % 3. alpha12; A52; ALGEBRAIC(:,27)
    alg(:,3) = 3.802./( 0.1027.*exp(-(alg(:,25)+x(3))./15.0)+ 0.23.*exp(-(alg(:,25)+x(3))./150.0));
    % 4. beta12; A55; ALGEBRAIC(:,38)
    alg(:,4) =  0.2.*exp(-(alg(:,25)-x(4))./20.3);
    % 5. alpha13; A53; ALGEBRAIC(:,32)
    alg(:,5) = 3.802./( 0.1027.*exp(-(alg(:,25)+x(5))./12.0)+ 0.25.*exp(-(alg(:,25)+x(5))./150.0));
    % 6. beta13; A56; ALGEBRAIC(:,40)
    alg(:,6) =  0.22.*exp(-(alg(:,25)-x(6))./20.3);
    % 7. alpha111
    alg(:,7) = 3.802./(0.1027.*exp(-(alg(:,25)+x(7))./17.0) + 0.2.*exp(-(alg(:,25)+x(7))./150.0));
    % 8. alpha112
    alg(:,8) = 3.802./(0.1027.*exp(-(alg(:,25)+x(8))./15.0) + 0.23.*exp(-(alg(:,25)+x(8))./150.0));
    % 9. beta111
    alg(:,9) = 0.1917.*exp(-(alg(:,25)+x(9))./20.3);
    % 10. beta112
    alg(:,10) = 0.2.*exp(-(alg(:,25)-x(10))./20.3);
    % 11. alpha31; A57; ALGEBRAIC(:,42)
    alg(:,11) = 7.0e-07.*exp(-(alg(:,25)+x(11))./x(23));
    % 12. beta31; A58; ALGEBRAIC(:,44)
    alg(:,12) = x(12)+ 2.0e-05.*(alg(:,25)+7.0);
    % 13. alpha32
    alg(:,13) = 7.0e-07.*exp(-(alg(:,25)+x(13))./x(24));
    % 14. beta32
    alg(:,14) = x(14) + 2.0e-05.*(alg(:,25)+7.0);
    % 15. alpha33
    alg(:,15) = 7.0e-07.*exp(-(alg(:,25)+x(15))./x(25));
    % 16. beta33
    alg(:,16) = x(16) + 2.0e-05.*(alg(:,25)+7.0);
    % 17. alpha2; A59; ALGEBRAIC(:,46)
    alg(:,17) = 1.0./( 0.188495.*exp(-(alg(:,25)+x(17))./x(22))+0.393956);
    % 18. beta2; A60; ALGEBRAIC(:,48) alpha13.*alpha2.*alpha33./(beta13.*beta33)
    alg(:,18) = (alg(:,5).*alg(:,17).*alg(:,15))./(alg(:,6).*alg(:,16));
    % 19. alpha4; A61; ALGEBRAIC(:,50)
    alg(:,19) = x(18).*alg(:,9);
    % 20. beta4; A62; ALGEBRAIC(:,52)
    alg(:,20) = x(19).*alg(:,11);
    % 21. alpha5; A63; ALGEBRAIC(:,54)
    alg(:,21) = x(20).*alg(:,9);
    % 22. beta5; A64; ALGEBRAIC(:,56)
    alg(:,22) = x(21).*alg(:,11);
    % 23. CNa3 (C3); A42; ALGEBRAIC(:,4)
    alg(:,23) = 1.0 - (states(:,1)+states(:,2)+states(:,3)+states(:,6)+states(:,4)+states(:,5)+states(:,7)+states(:,8));
    % 24. I_Na; A40; ALGEBRAIC(:,58)
    alg(:,24) =  constants(1).*states(:,1).*(alg(:,25) - constants(2));
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
