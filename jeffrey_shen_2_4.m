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


H = diag([En En]);
V = 50e-3 * 1.602677e-19; % Joules
rho =   sqrt(2*m*(V-En)/(hbar^2));
T = 16*En.*(V-En).*exp(-rho*l)/V^2;
delta = hbar *(sqrt(2*En/m)/(2*a)).*T;
toeV = 6.241509e18;
H = H + diag(delta,2) + diag(delta,-2);


len_L = length(x_L);
len_R = length(x_R);
len_X = len_L + len_R;
x = [x_L x_R];

b1 = [phi_1L zeros(1,len_R)];
b2 = [phi_2L zeros(1,len_R)];
b3 = [zeros(1,len_L) phi_1R];
b4 = [zeros(1,len_L) phi_2R];

psi_0 = [1 1 0 0]' * (1/sqrt(2));


time = 0:1e-15:1e-11;
T = length(time);

psi_t = zeros(T,4);

psi_t_x = zeros(T,len_X);


for i = 1:T
   U = expm(-1i * H * time(i)/hbar);
   psi_t(i,:) = U*psi_0;
   psi_t_x(i,:) = psi_t(i,1)*b1 + psi_t(i,2)*b2;
end


figure(1)
imagesc(x,time,abs(psi_t_x).^2);
shading flat
colorbar;
xlabel('position');
ylabel('time');
title('Probability Density |\Psi_t|^2');

