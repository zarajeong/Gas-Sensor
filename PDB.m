% Persamaan Differensial

function [dxdt] = PDB(tt,x,T,PH2S,PO,Ef2)

kT = 8.625e-5*(165+273); %Boltzmann
sigmaPt = 1.3E+19; %material Pt
mH2S = 5.66e-26; %massa H2S
mO = 2.6560e-26; %massa O
mSO2 = 1.0240e-25; %massa SO2
mH2O = 2.98897e-26; %massa H2O

% TotalX adalah total coverage
TotalX = x(1) + x(2) + x(3) + x(4); %H2S + 3O -> SO2 + H2O 
% Material yang digunakan sementara adalah Pt (k = 1.38E-23 J)
fH2S = (1-TotalX)*PH2S*1.01E+5/(sigmaPt*sqrt(2*pi*mH2S*1.38E-23*T)); %(2.12)
fO = ((1-TotalX)^2)*PO*1.01E+5/(sigmaPt*sqrt(2*pi*mO*1.38E-23*T)); %(2.14)

% Untuk rf1 (H2S)
Sx1 = 2.28;  
Ef1 = 0.084;    
rf1 = Sx1*exp(-Ef1/kT);

% Untuk rr1 (H2S)
vr1 = 9.5e+12;  
Er1 = 0.71;     
rr1 = vr1*exp(-Er1/kT);

% Untuk rf2 (O)
Sx2 = 0.51;  
% Ef2 = 0.015;  
rf2 = Sx2*exp(-Ef2/kT);

% Untuk rr2 (O)
vr2 = 1.4e+12;    
Er2 = 0.174;     
rr2 = vr2*exp(-Er2/kT);

% Untuk rf3 (H2S dan O)
vf3 = 1.4e+12;     
Ef3 = 0.284;   
rf3 = vf3*exp(-Ef3/kT);

% Untuk rr3 (H2S dan O)
vr3 = 1.7e+12;      
Er3 = 0.086; 
rr3 = vr3*exp(-Er3/kT);

% Nilai differensial
dH2S_dt = rf1*fH2S - rr1*x(1) + rr3*x(3)*x(4) - rf3*x(1)*x(2)^3;
dO_dt   = 2*rf2*fO - 2*rr2*x(2)^2 + rr3*x(3)*x(4) - rf3*x(1)*x(2)^3;
dSO2_dt = rf3*x(1)*x(2)^3 - rr3*x(3)*x(4);
dH2O_dt = rf3*x(1)*x(2)^3 - rr3*x(3)*x(4);
dxdt    = [dH2S_dt; dO_dt; dSO2_dt; dH2O_dt];
end