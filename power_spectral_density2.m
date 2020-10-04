% Plot of Power Spectral Density of initial conditions

clear workspace;
close all;

% First Initial Condition
figure;
dx = 0.01;
Fs = 1/dx;  % Sampling Frequency
x = 0:dx:0.5-dx;
y = (sin(2*pi*x)).^2;
[Pxx,F] = periodogram(y,[],length(y),Fs,'onesided'); 
stem(F,Pxx, 'linewidth', 1.5);
xlabel('Normalizd Frequency, \beta\Deltax');
ylabel('power');
title('Power Spectral Density, FFT length (N) = 100');

% Second Initial Condition
figure;
dx = 0.01; 
Fs = 1/dx;  % Sampling Frequency
x = 0:dx:0.5-dx;  
y = (sin(2*pi*x)).^20;
[Pxx,F] = periodogram(y,[],length(y),Fs,'onesided');
stem(F/length(y),Pxx, 'linewidth', 1.5);
xlabel('Normalizd Frequency, \beta\Deltax');
ylabel('power');
title('Power Spectral Density, FFT length (N) = 100');


