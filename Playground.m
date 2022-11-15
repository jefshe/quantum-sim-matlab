clear all
clc

% this code calculates the time evolution of an ammonia molecule, truncated
% to its two lowest energy states. We use as basis states |phi_L>, |phi_R>,
% i.e. the localized modes in the left and right potential well.
% A tunneling matrix element Delta describes the transparency of the
% barrier
% A detuning epsilon describes the effect of an electric field that pulls
% the molecule towards one side


%% Definitions and parameters

hbar = 1; 
% I am setting the Planck constant to 1, which means I am expressing the
% energies in frequency units

% Pauli matrices

sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];  % 1i is the imaginary unit
sigma_z = [1 0; 0 -1];   

% Basis states. These are the 2-dimensional representations of the "real
% state". In other words, they are the projections in the 2x2 subspace

pL = [1; 0];
pR = [0; 1];

% Parameters of the potential well (see lecture notes for sketch)
a = 3;   % width of each well
b = 2;  % displacement of the bottom of each well

% Tunneling and detuning

Delta = 4*pi*2e6;   % tunnel matrix element
epsilon = 4*pi*0e6; % detuning

% Hamiltonian

H = Delta * sigma_x + epsilon * sigma_z;

% range of x-coordinates
x = -4:0.01:4;
L = length(x);

%% Basis states

phi_L = zeros(1,L);
phi_R = zeros(1,L);

for i=1:L
    if (x(i) > (-b - a/2)) && (x(i) < (-b + a/2))
        phi_L(i) = sqrt(2/a) * sin(pi * (x(i) - (-b-a/2)) / a);
    end
    if (x(i) > (b - a/2)) && (x(i) < (b + a/2))
        phi_R(i) = sqrt(2/a) * sin(pi * (x(i) - (b-a/2)) / a);
    end
end

figure(1)
plot(x,abs(phi_L).^2,x,abs(phi_R).^2);  
% plot the square of the wavefunctions --> probability of finding the particle in some point

% MATLAB tip: when you add a "." before an operation (the square in this
% case) it does the operation element by element. If you just wrote "^2"
% you would get an error, because Matlab would then try to do the square of
% matrix, but this array is not a square matrix

%% Time evolution

time = 0:1e-8:1e-6;
T = length(time);

% Preallocate the array that will contain on the columns the wavefunction at
% time t. This array has L columns, as many as the x-points we used to
% describe the wavefunctions, and T rows, as many points in time we are
% going to evaluate

psi_t = zeros(T,L);

% preallocate an array of the coefficients of the basis states, which will
% be used to express the actual wavefunctions. Note: we defined the
% coefficient vectors as column vectors, so this array must have 2 columns
% and T rows

p_t = zeros(2,T);

% Initial state = |phi_L>

p0 = pL;
% p0 = 1/sqrt(2)*(pL+pR);

for i = 1:T
   U = expm(-1i * H * time(i));   
   % time evolution operator. Notice that there is not hbar, since I set it
   % to 1 from the beginning
      
   p_t(:,i) = U * p0; % time-evolved coefficients
   psi_t(i,:) = phi_L * p_t(1,i) + phi_R * p_t(2,i);
end

% plot the results. I am using a color-scale plot to show the square of teh
% wavefunction as a function of both position and time.
% MATLB tip: the command "imagesc" does a color-plot of a function of two
% variables
figure(2)
imagesc(x,time,abs(psi_t).^2);
shading flat
colorbar;
xlabel('position');
ylabel('time');
title('Probability Density |\Psi_t|^2');


%% Exercise for next week:
% Familiarize yourself with this code. Play with the parameters and learn
% how to optimize the looks of the output figure. Matlab gives many more
% options to plot data of this type! 



