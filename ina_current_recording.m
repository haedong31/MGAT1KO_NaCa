clc
clearvars
close all

% initial parameter values
x0 = [2.5,2.5,2.5,-2.5,2.5,-7.5,2.5,2.5,2.5,-2.5,7,0.0084,7,0.0084,7,0.0084,...
    7,1/1000.0,1.0,1/95000.0,1/50.0,16.6,7.7,7.7,7.7,13.0];

% experimental protocol
prot = NaN(7,1);
prot(1) = 39.4843; % ENa
prot(2) = -100; % holding potential
prot(3) = 5; % holding time
prot(4) = 5; % P1
prot(5) = 120; % P1 time
prot(6) = -100; % P2
prot(7) = 25; % P2 time

[t,s,a] = INa(x0,prot);

figure('Color','w')
plot(t,a(:,24),'Color','red','LineWidth',1.5)
title("I_{Na}")
