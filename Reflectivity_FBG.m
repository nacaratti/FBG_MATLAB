%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              DAVI PONTES NACARATTI - 03/05/2021 (Last review)          %
%   Modeling of Fiber Bragg Gratings with Different Lengths for the      %
%                  Reflectivity Control for Fiber Lasers                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

L = 25e-3;                % comprimento da rede
neff = 1.45;              % indice de refracao efetivo
p = 366.88e-9;            % periodo da rede de bragg
dn = 6e-5;                % variacao do indice de refracao
lambdab = 2*neff*p;       % comprimento de onda (condição de Bragg)
v = 0.7;                  % visibilidade das franjas

lambda1 = lambdab-5e-10:1e-12:lambdab+5e-10;
lambdaN = (lambdab - lambda1)/lambdab;
lmax = (1+dn/neff)*lambdab;

for ii=1:length(lambda1)
    lambda = lambda1(ii);
    k = pi*v*dn/lambda;
    db = 2*pi*neff*(1/lambda - 1/lambdab);  
    sigma = 2*pi*dn/lambda;
    sigma1 = db + sigma;
    gamma = sqrt (k^2 - sigma1^2);
    DL(ii) = (lambda^2/(pi*neff*L))*sqrt(k^2*L^2+pi^2); % Largura de banda: Eq. 4.6.14 Raman Kashyap
    
    % Matriz de Transferencia
    T11 = cosh(gamma*L) - 1i*sigma1*sinh(gamma*L)/gamma;
    T22 = cosh(gamma*L) + 1i*sigma1*sinh(gamma*L)/gamma;
    T12 = -1i*k*sinh(gamma*L)/gamma;
    T21 =  1i*k*sinh(gamma*L)/gamma;
    
    R(ii) = abs(T21/T11)^2; % Refletividade
    S(ii) = abs(1/T11)^2;   % Transmissividade
end

Rm = max(R);
Dl = max(DL);

plot(lambda1,R);
xlabel('Wavelength(m)');
ylabel('Reflectivity (x 100%)');
title('');
grid on
hold on