% 1kT = phi = 0.125; Lpi = 0.0077
% 2kT = phi = 0.202; Lpi = 0.0149
% 3kT = phi = 0.263; Lpi = 0.0220
% 4kT = phi = 0.314; Lpi = 0.0289
% 5kT = phi = 0.360; Lpi = 0.0359
clear
clc
close all

phi = 0.36;
a = 1435;
T = 20;
L = a/5;
v = 4*pi*(L*1e-9)^3/3;
n = phi/v; % m^-3
kT = 1.38*10^-23 * (273+T); %kg*m^2*s^-2
OSMOTIC = n*(1+4*phi); %m^-3

pfpp = 2.2975*a;
e = 1.6e-19;
eps = 8.854e-12*80;
r = 0+2*a:2*a+2*L;
kappa = 143.5;

B = 32*pi*eps*(a*1e-9)*(kT/e)^2*(tanh(e*50*10^-3/4/kT))^2; %SI Unit
BC = 1e18*kT*pfpp;
RepPot = B*exp(-(r-2*a)/10); %SI
RepForceSI = B*10^8*exp(-(r-2*a)/10);
RepNor = RepPot/kT;

LPi = OSMOTIC*kT; %SI Unit
% LPi = 0.0612;
AttPot = -(4/3)*LPi*pi*((a+L)*1e-9).^3*(1-0.75*r/(a+L)+0.0625*(r).^3/((a+L).^3));
AttForceSI = -(4/3)*LPi*pi*(-0.75*((a+L)*1e-9)^2+3*(r.*1e-9).^2/16);
AttFOrceC = AttForceSI*1e9;
AttNor = AttPot/kT;

Potential = RepNor + AttNor;
Pmin = min(Potential)

% figure(1)
% hold on
% plot(r-2*a,RepNor,'b--')
% plot(r-2*a,AttNor,'r--')
% plot(r-2*a,Potential,'k')
% plot(0:2*L,zeros(1,2*L+1),'k--')

% axis([0, 2*L  -10 1])
% xlabel('r (nm)')
% ylabel('U(r)/kT')
% legend('Repulsive','Attractive','Total Potential')