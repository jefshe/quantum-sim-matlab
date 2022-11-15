clear all
clc
close all

%% Definitions and parameters
hbar = 6.62606876e-34 / (2*pi); %Planck constant
m = 9.10938188e-31;  % electron mass
a = 6e-9;  % width of the well
l = 2e-9;

n = 1:2;
En = (n*pi*hbar).^2 / (2*m*a^2);

dx = 1e-11;
x_L = -l-a:dx:-l;
x_R = l:dx:l+a;

phi_1L = sqrt(2/a) * sin(pi * (x_L + l+a)/ a);
phi_2L = sqrt(2/a) * sin(2*pi * (x_L + l+a) / a);

phi_1R = sqrt(2/a) * sin(pi * (x_R - l) / a);
phi_2R = sqrt(2/a) * sin(2*pi * (x_R - l) / a);



%% Plot the results.    
figure(1)
%subplot(2,2,1)
plot(x_L,abs(phi_1L).^2,x_L,abs(phi_2L).^2,x_R,abs(phi_1R).^2,x_R,abs(phi_2R).^2); 
title('Stationary states of 2 infinite potential wells');
xlabel('x (nm)');
ylabel('P(x)');
legend('phi 1L','phi 2L','phi 1R','phi 2R')