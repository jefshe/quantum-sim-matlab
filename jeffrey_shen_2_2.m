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
disp('E1 in eV:')
disp(num2str(En(1) * 6.241509e18))
disp(' ')
disp('E2 in eV:')
disp(num2str(En(2) * 6.241509e18))
