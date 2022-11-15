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
H = diag([En En]);
V = 50e-3 * 1.602677e-19; % Joules
rho =   sqrt(2*m*(V-En)/(hbar^2));
T = 16*En.*(V-En).*exp(-rho*l)/V^2;
delta = hbar *(sqrt(2*En/m)/(2*a)).*T;
toeV = 6.241509e18;

H = H + diag(delta,2) + diag(delta,-2);

disp('H in meV:')
disp(' ')
disp(num2str(H*toeV*1000))
disp(' ')
disp('With eigen energies of:')
disp(' ')

[vec,val] = eig(H*toeV*1000);

disp(eig(H*toeV*1000))
disp(' ')
disp('With corresponding eigen vectors of:')
disp(' ')
disp(num2str(vec))

