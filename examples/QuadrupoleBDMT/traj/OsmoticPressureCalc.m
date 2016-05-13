clear
clc
close all
i =0;
% for phi = 0.1:0.1:0.5
    phi = 0.4;
    a = 1435;
    T = 20;
    L = a/5;
    v = 4*pi*(L*1e-9)^3/3;
    n = phi/v; % m^-3
    kT = 1.38*10^-23 * (273+T); %kg*m^2*s^-2
    OSMOTIC = n*(1+4*phi); %m^-3
    OSnm = OSMOTIC*10^-27;
    
    pfpp = 2.2975*a;
    e = 1.6e-19;
    eps = 8.854e-12*80;
    r = 0+2*a:2*a+2*L;
%     TT = 0+2*a:3*a+2*a;
%     r = r+TT
%     AttNor = 0*TT;
%     r = 4305.4304999999986;
    kappa = 143.5;
    
    B = 32*pi*eps*(a*1e-9)*(kT/e)^2*(tanh(e*50*10^-3/4/kT))^2; %SI Unit
    LPi = OSMOTIC*kT; %SI Unit
%     LPi = 0.6355;
    BC = 1e18*kT*pfpp;
    RepPot = B*exp(-(r-2*a)/10); %SI
    RepForceSI = B*10^8*exp(-(r-2*a)/10);
    RepForceC = RepForceSI*1e9;
    RepNor = RepPot/kT;
    hold on
    plot(r-2*a,RepNor,'b--')
    % plot(Rep)
    AttPot = -(4/3)*LPi*pi*((a+L)*1e-9).^3*(1-0.75*r/(a+L)+0.0625*(r).^3/((a+L).^3));
    AttForceSI = -(4/3)*LPi*pi*(-0.75*((a+L)*1e-9)^2+3*(r.*1e-9).^2/16);
    AttFOrceC = AttForceSI*1e9;
    AttNor = AttPot/kT;
%     AttNor = [AttNor, zeros(1,3511-188)];
    plot(r-2*a,AttNor,'r--')
    Potential = RepNor + AttNor;
    % figure(2)
    plot(r-2*a,Potential,'k')
    axis([0, 2*L  -10 1])
    xlabel('nm');
    [aa bb] = min(Potential);
    title(num2str(bb));

% end
%     xlabel('r (nm)')
%     ylabel('U(r)/kT')
%     legend('Repulsive','Attractive','Total Potential');
%     legend('\phi = 0.1','\phi = 0.2','\phi = 0.3','\phi = 0.4','\phi = 0.5')
    plot(0:200,zeros(1,201),'k--')
% 

