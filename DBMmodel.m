function Absorption = DBMmodel(sigma,thick,freq)

rho0 =1.204; %Denotes value in air where ambiguity might otherwise arise [kg/m3]
c0 =343; %Speed of sound (m/s)
z0 =rho0*c0; %Specific acoustic impedance of air (Kg/m2*s)
omega =2*pi*freq; %Angular frequency [s^-1]


z_Miki = rho0*c0*( 1 + 5.50*((freq./sigma)*1000).^(-0.632)- 1i*8.43*((freq./sigma)*1000).^(-0.632) );
k_Miki = omega/c0 .* (-i) .* ( 11.41*((freq./sigma)*1000).^(-0.618)+ 1i* (1 + 7.81*((freq./sigma)*1000).^(-0.618) ) );
%% Compute sound absorption for a sample of thickness d backed by a rigid and impervious wall under at room temperature and pressure conditions
Zs = -i.*z_Miki./tan(k_Miki*thick);

Absorption = 1 - ( abs( (Zs-rho0*c0)./(Zs+rho0*c0) ) ).^2;