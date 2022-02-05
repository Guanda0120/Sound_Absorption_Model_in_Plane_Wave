function Absorption = JCAmodel(tort,sigma,epsilon,viscous,thermal,thick,freq)

%% Varible Parameter Setting
%tort:      Tortuosity [-] range[1-3]
%sigma:     Airflow resistivity [Ns/m4]
%epsilon:   Porosity [-] range [0-1]
%viscous:   Viscous characteristic dimension [m] [range 10-1000 10^-6m]
%thermal:   Thermal characteristic dimension [m] [range 10-1000 10^-6m]
%thick:     thickness [m]
%freq:      frequency

%% Constant Parameter Setting
rho0 =1.204;        %Denotes value in air where ambiguity might otherwise arise [kg/m3]
c0 = 343;           %Speed of sound (m/s)
z0 =rho0*c0;        %Specific acoustic impedance of air (kg/m2*s)
eta = 1.84*10^-5;   %Visocity of air [Poiseuille] (1.84*10^-5 )
omega =2*pi*freq;   %Angular frequency [s^-1]
gamma =1.4;         %Ratio of the specific heat capacity [-]
Npr =0.71;          %Prandtl number (0.77 at 20*C)
P0 =101320;         %Atmospheric pressure [N/m2]

%% Computation
%Accounts for the viscous losses
rho_eq=(tort*rho0/epsilon).*(1+((sigma*epsilon)./(1i*omega*rho0*tort)).*sqrt(1+((1i*4*(tort^2)*eta*rho0*omega)/((sigma^2)*(viscous^2)*(epsilon^2)))));

%Accounts for the thermal losses
K_eq=(gamma*P0/epsilon)./(gamma-(gamma-1).*(1+(((8*eta)./(((thermal*10^-6)^2)*Npr*rho0*omega*1i)).*sqrt(1+(1i*(((thermal*10^-6)^2)*Npr*rho0*omega)/(16*eta))))).^-1);%Characteristic Impedance

Zc=(K_eq.*rho_eq).^0.5;
%Complex wave number
k_c=omega.*sqrt(rho_eq./K_eq);
%Surface impedance of sample
Zs=-1i*Zc.*cot(k_c*thick);
%Normalised specific acoustic impedance
Zn = Zs./z0;
%Reflection coefficient
R=(Zs-z0)./(Zs+z0);
%Absorption coefficent
Absorption = 1-(abs(R)).^2;
end
