function Absorption = AttenModel(tort,sigma,epsilon,d,fvector)
%d =0.05; %Thickness of the sample

%fvector =2:2:6400; %Frequency range
rho_0 =1.204; %Denotes value in air where ambiguity might otherwise arise [kg/m3]
c0 =343; %Speed of sound (m/s)
z0 =rho_0*c0; %Specific acoustic impedance of air (Kg/m2*s)
%tort =1.5; %Tortuosity [-] range[1-3]
%sigma =62900; %Airflow resistivity [Ns/m4]
%epsilon =0.28; %Porosity [-]
omega =2*pi*fvector; %Angular frequency [s^-1]
eta =18.27e-6; %Visocity of air [Poiseuille] 1.84*10^-5
P0 =101320; %Atmospheric pressure [N/m2]
gamma =1.4; %Ratio of the specific heat capacity [-]
Cp =1.01; %Specific heat capacity of air at constant pressure [J/kg/K]
kappa =2.41*10^-2; %Thermal conductivity of air [Wm/K]
Npr = 0.71; %Prandtl number (0.77 at 20*C)
s_b = 1.0; %Shape factor
% Effective (/dynamic) Density
% Accounts for the viscous losses
%-----------------------------------------------------------------

lambda=(1/2*s_b)*((8*tort*rho_0*omega)./(sigma*epsilon)).^0.5;
T_bes=besselj(1,lambda*sqrt(-i))./besselj(0,lambda*sqrt(-i));
rho_c=((tort*rho_0)/epsilon)*(1-(2./lambda*sqrt(-i)).*T_bes).^-1;

% Effective (/dynamic) bulk modulus
% Accounts for the thermal losses
%-----------------------------------------------------------------
T_bes2=besselj(1,Npr^0.5*lambda*sqrt(-i))./besselj(0,Npr^0.5*lambda*sqrt(-i));
K_w=((gamma*P0)/epsilon)*(1+((2*(gamma-1))./(Npr^0.5*lambda*sqrt(-i))).*T_bes2).^-1;
%Characteristic Impedance
Zc=(K_w.*rho_c).^0.5;
%Complex wave number
k_c=omega.*sqrt(rho_c./K_w);
%Surface impedance of sample
Zs=-1i*Zc.*cot(k_c*d);
%Normalised impedance
Zn=Zs/z0;
%Reflection coefficient
R=(Zs-z0)./(Zs+z0);
%Absorption coefficent
Absorption=1-(abs(R)).^2;
