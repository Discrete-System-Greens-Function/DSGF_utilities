%function [ epsilon ] = SiC_dielectric_function( omega )

% This code uses the Drude-Lorentz model to calculate the dielectric
% function of silicon carbide (SiC).

% INPUTS:  omega         vector of frequencies at which to calculate the dielectric function [rad/s]
%
%
% OUTPUTS: epsilon       dielectric function [-]
back = cd;
saveDir = fullfile(back,'../../library/Non_uniform_spectra/');

%omega = readmatrix([saveDir 'SiC_non_uniform_spectra_85.csv']);
omega = readmatrix(['SiC_non_uniform_spectra_85.csv']);
% Dielectric function approximated with a Lorentz model
omega_TO = 1.494e14;              % Transverse optical phonon frequency [rad/s]
omega_LO = 1.825e14;              % Longitudinal optical phonon frequency [rad/s]
Gamma_D = 8.966e11;               % Damping constant [rad/s]
epsilon_inf = 6.7;                % High-frequency limit of permittivity
epsilon = epsilon_inf*( ((omega.^2) - (omega_LO^2) + 1i*Gamma_D*omega)./ ((omega.^2) - (omega_TO^2) + 1i*Gamma_D*omega) );

% Complex refractive indices
%m = sqrt(epsilon);

epsilon_fig = figure(1);
semilogx(omega, real(epsilon), omega, imag(epsilon), '--', 'linewidth', 2) % '--'
%plot(lambda*1e6, real(epsilon), lambda*1e6, imag(epsilon), '--', 'linewidth', 2)
xlabel('Frequency [rad/s]', 'fontsize', 12)
%xlabel('Wavelength, \lambda [\mum]', 'fontsize', 12)
ylabel('Dielectric function, \epsilon', 'fontsize', 12)
title(['Dielectric function of ' material ', N_o_m_e_g_a = ' num2str(N_omega)], 'fontsize', 16)
legend('Real part', 'Imaginary part', 'location', 'best')
set(gca, 'fontsize', 16)
axis tight
grid on

saveas(epsilon_fig, 'fig_SiC_dieletric_function.fig');
saveas(epsilon_fig, 'fig_SiC_dieletric_function.png');