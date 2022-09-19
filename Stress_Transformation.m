%%
%%% Arian Velayati, PhD
%%%% This script is used to find stresses in new coordinates 
... using the stress transformation tensor method
... The procedure to calculate stresses (sigma_n,tau_n)on an arbitrary plane
... given its orientation respect to the geographical coordinate system (strike,dip)
... and the in-situ stress tensor of principal stresses Sp (given its principal values and principal directions).
    clc; close; clear;
%% Input, principal stresses %%%% ALL stressa units can be both MPa and psi
Sp =[2400 0 0;0 1200 0;0 0 1000]; %MPa
Shmin_azi = 90;
SHmax_azi = 150;
Pp = 440; %MPa
strike = 120; dip = 70;

%% Director angles
I = input('Select the faulting regime: NF = 1; SS = 2; RF = 3:  ');

if I == 1
% Normal Faulting
a = Shmin_azi;
B = 90;
Gamma = 0;
    elseif I == 2
a = SHmax_azi;
B = 0;
Gamma = 90;
    else
a = SHmax_azi;
B = 0;
Gamma = 0; 
end

RPG = [cosd(a)*cosd(B), sind(a)*cosd(B), -sind(B);
       cosd(a)*sind(B)*sind(Gamma)-sind(a)*cosd(Gamma), sind(a)*sind(B)*sind(Gamma)+cosd(a)*cosd(Gamma), cosd(B)*sind(Gamma);
       cosd(a)*sind(B)*cosd(Gamma) + sind(a)*sind(Gamma), sind(a)*sind(B)*cosd(Gamma)-cosd(a)*sind(Gamma), cosd(B)*cosd(Gamma)];

SG = RPG' * Sp * RPG
%% The coordinate system basis is comprised of nd (dip), ns (strike), and nn (normal) vectors:
... d-s-n right-handed basis. The three vectors depend solely in two variables: strike and dip of the fault.
    

% Fault coordinate system as a function of strike and dip
nn = [-sind(strike)*sind(dip); cosd(strike)*sind(dip); -cosd(dip)];
ns = [cosd(strike); sind(strike); 0];
nd = [-sind(strike)*cosd(dip); cosd(strike)*cosd(dip); sind(dip)];

% The fourth step consists in projecting the stress tensor based on the geographical coordinate system onto the fault base vectors. 
... The stress vector acting on the plane of the fault is t ({t} is not necessarily aligned with nd, ns, or nn):

t = SG*nn %MPa

% The total normal stress on the plane of the fault is Sn (aligned with nn)

Sn = dot(t,nn) %MPa

% The effective normal stress on the fault plane is:

sig_n = Sn - Pp;

% The shear stress on the plane of fault is aligned with nd and ns are:

tau_d = dot(t,nd)
tau_s = dot(t,ns)

% Effective norma stress and absolute shear can also be calculated with the
% following eqns

Sig_N = dot(t,nn)-Pp
Tau_abs = sqrt(tau_d^2 + tau_s^2)


% The rake is the angle of the shear stress tau_d+tau_s with respect to ns (horizontal line) and quantifies the direction of expected fault movement in the fault plane.
    
rake = atand(tau_d/tau_s)
   