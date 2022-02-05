function Absorption = H&S_model(tort,sigma,phy,deviationPoreSize,h,freq)
% tort = 1.45;
% sigma = 3663;
% phy = 0.45;
% deviationPoreSize = 0.74;
% h = 0.0965;
% freq = 1:1:1600;

rho_0 =1.204;        %Denotes value in air where ambiguity might otherwise arise [kg/m3]
c0 = 343;           %Speed of sound (m/s)
garmma =1.4;         %Ratio of the specific heat capacity [-]
Npr =0.71;          %Prandtl number (0.77 at 20*C)
P0 =101320;         %Atmospheric pressure [N/m2]
omega = 2*pi.*freq;

delta = (deviationPoreSize*log(2))^2;
theta_1 = (4/3)*exp(4*delta)-1;
theta_2 = (1/sqrt(2))*exp(1.5*delta);
a_1 = theta_1/theta_2;
a_2 = theta_1;
b_1 = a_1;
epsilon = sqrt(-i.*omega*rho_0*(tort^2)/(sigma*phy));
Fomega = (1+a_1.*epsilon+a_2*exp(2))./(1+b_1.*epsilon);
Fomega_Npr = (1+a_1.*epsilon*Npr+a_2*exp(2))./(1+b_1.*epsilon*Npr);

rho_x = 1+Fomega./((tort.*epsilon).^2);
rho_b = (rho_0*tort^2/phy).*rho_x;
rho_xNpr = 1+Fomega_Npr./((tort.*epsilon*Npr).^2); 
ci_b = phy*(1/(garmma*P0))*(garmma-(garmma-1)./rho_xNpr);

Z_b = sqrt(rho_b./ci_b)./(rho_0*c0);
K_b = omega.*sqrt(rho_b.*ci_b);

Z_s = Z_b.*coth(-i.*K_b*h);
Absorption = 1-abs((Z_s-1)./(Z_s+1)).^2;
