%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System ODE Script using dsolve function to solve FBG Reflectivity and
% Transmissivity
% Author: Davi Pontes Nacaratti
% Date: 20/01/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

syms R(z) S(z) s(lambda) k(lambda) delta(lambda) lambda r(lambda)

L = 10e-3;          % FBG Length
neff = 1.45;        % Core Refractive index
dneff = 1e-4;       % Refractive index variation
lambdaB = 1018e-9;  % Center wavelength
v = 1;              % Visibility (0<v<1)

delta(lambda) = 2*pi*neff*(1/lambda - 1/lambdaB);
k(lambda) = pi*v*dneff/lambda;
sigma(lambda) = 2*pi*dneff/lambda;
s(lambda) = delta + sigma;

% System of ODE eqs. 
eqns = [diff(R,z) == 1i*s(lambda)*R + 1i*k(lambda)*S; diff(S,z) == -1i*s(lambda)*S - 1i*k(lambda)*R];
% Boundary conditions
cond = [R(0)==1 ; S(L)==0];
% Solve R and S
[RSol(z),SSol(z)] = dsolve(eqns,cond);
% Reflectivity
r(lambda) = SSol(0)/RSol(0);
lambda = linspace(1016e-9,1020e-9,1000);
solr = abs(double(r(lambda)));
solt = 1-abs(double(r(lambda)));
%solr = 10*log10(abs(double(r(lambda))));
%solt = 10*log10(1-abs(double(r(lambda))));
% Ploting
plot((lambda-lambdaB)/10e-9,solr);
hold on;
plot((lambda-lambdaB)/10e-9,solt);
xlabel('Comprimento de onda Normalizado (nm)');
ylabel('Refletividade (%)');
%ylabel('Refletividade (dB)');
title('Gráfico Refletividade');
grid on