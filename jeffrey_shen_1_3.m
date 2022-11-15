clear all
clc
close all

%% Definitions and parameters
hbar = 6.62606876e-34 / (2*pi); %Planck constant
S_x = (1/2) * [0 1; 1 0];
S_y = (1/2) * [0 -1i; 1i 0];  % 1i is the imaginary unit
S_z = (1/2) * [1 0; 0 -1];      
plus = [1; 0];
minus = [0; 1];
B_0 = 1;     % Static magnetic field = 1T
B_1 = 2e-4;  % Oscillating magnetic field, at frequency omega_0
gamma = 2*pi*28e9;  % electron gyromagnetic ratio = 28 GHz/T
omega_0 = gamma * (B_0);   % Larmor frequency of a spin in the field B_0
omega_1 = gamma * B_1;   % Rabi frequency. This is the frequency at which a 
                         % spin would precess if placed in a static field B_1
figure(1)  
i = 1;
%for noisePower = [0.01e-6 0.1e-6 1e-6 10e-6 100e-6]                       
for noisePower = [0.001 0.005 0.01 0.05]                       
    %% Time evolution
    dt = 3e-12;  % time sampling interval. How do you choose this value properly?
    tau = 0:3e-12:3.57e-11/2;
    avgSz=zeros(length(tau),1);
    H0 = (omega_0) * S_z;

    half_period = 8.93e-8;
    psi0 = minus;   %initial state.
    time = 0:dt:half_period;

    B_n = noisePower*wgn(1,length(time),0);
    for t = 1:length(time)
       % Hamiltonian describing the rotating magnetic field
       %H1 = omega_1 * (S_x * cos(omega_0 * time(t)) + S_y * sin(omega_0 * time(t))); 
       H1 = omega_1 * (S_x * cos(omega_0 * time(t))); 
       U = expm(-1i * (H0 + H1 + B_n(t)*gamma*S_z) * dt);   % time evolution operator. No need to divide 
       % by hbar, since I wrote the Hamiltonian in frequency units. The
       % Hamiltonian in there is the sum of the static and the rotating parts
       psi0 = U * psi0; % time-evolved state
    end

    origPsi = psi0;
    time2 = half_period:dt:2*half_period;

    for u  = 1:length(tau)
        psi0 = origPsi;

        dt2 = 3e-13;
        rotateTime = 0:dt2:tau(u);

        %Free Evolution
        B_n = noisePower*wgn(1,length(rotateTime),0);

        for t = 1:length(rotateTime)
            U = expm(-1i * (H0 + B_n(t)*gamma*S_z) *dt2);
            psi0 = U*psi0;
        end
        %Rotate 90 Degrees
        B_n = noisePower*wgn(1,length(time2),0);
        
        for t = 1:length(time2)
               H1 = omega_1 * (S_x * cos(omega_0 * time2(t))); 
               U = expm(-1i * (H0 + H1 + B_n(t)*gamma*S_z) * dt);   % time evolution operator. No need to divide 
               % by hbar, since I wrote the Hamiltonian in frequency units. The
               % Hamiltonian in there is the sum of the static and the rotating parts
               psi0 = U * psi0; % time-evolved state
        end
        avgSz(u) = psi0' * S_z * psi0;

    end
    subplot(5,2,i)
    i= i+1;
    plot(tau, avgSz);
    ylim([-0.5, 0.5])
    xlabel ('Tau (s)');
    title(num2str(noisePower))
    ylabel('<Sz>');
    subplot(5,2,i)
    i = i+1;
    plot(time2,B_n)
end

%% Plot the results
% First panel: spin expectation value along X, Y, Z
%plot(time,avgSx,time,avgSy,time,avgSz)




