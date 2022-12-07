function [Time,Pressure,PressureOffset2,Temperature,Porosity,M_D,Dehyd,...
    DELTAG,Depth,P_tot,P_eff] = Dehydration(m_d0,Vars_in,T_in,T_end)
%%% This code solves the dimensional equations for pore pressure,
%%% temperature, porosity, and chemically bound water content, in a
%%% dehydrating and subducting 1-D column of material. For numerical 
%%% stability, an approximately 100 m long section of non-dehydrating
%%% material is included at the base of the column.

%%% This code employs Dirichlet conditions for pressure and temperature at 
%%% the bottom boundary. At the top boundary, it employs either a Dirichlet
%%% condition for temperature and an absorbing boundary condition for
%%% pressure/porosity.

%%% REQUIRED INPUT
% m_d0                        initial bound water content [kg m^(-3)]

%%% OPTIONAL INPUT FOR RESTARTING PREVIOUS SIMULATIONS
% Vars_in                     A 4*N by 1 column vector containing the
%                             output of a previous simulation. Must be
%                             defined as: 
%                             Vars_in = [pressure;temperature;phi;m_d];
% T_in                        Start time of continued simulation [s]
% T_end                       End time for continued simulation [s]

%%% OUTPUT
% Time                         Time vector of simulation [yr]
%%% The following output variables are arrays that give values at every
%%% node in the column for each time contained in the Time vector.
% Pressure                     Excess fluid pressure in the column [Pa]
% PressureOffset2              Fluid Pressure above the lithostatic fluid
%                              pressure [Pa]
% P_tot                        Total fluid pressure [Pa]
% P_eff                        Effective pressure [Pa] 
% Temperature                  Temperature in the column [K]
% Porosity                     Porosity in the column [-] 
% M_D                          Chemically released water content [kg m^(-3)]
% Dehydration                  Normalized dehydration rate [s^(-1)]
% DELTAG                       Gibbs free energy of the reaction [cal/mol]
% Depth                        Depth referenced to the top of the column [m]

%%% Mechanical parameters. Constants.
C_s = 2.5e6;        % volumetric heat capacity of rock [J m^(-3) K^(-1)]
k_i = 1e-20;        % background permeability [m^2] 1e-20
K_f = 2.9e9;        % fluid bulk modulus [Pa]
K_T = 2.25;         % thermal conductivity [W/(m K)
n = 3;              % perm-porosity power law exponent [-]
%%% P-wave velocity estimated form velocity model of Trehu et al. 2008.
%%% Audet et al. 2009 use a velocity model with mantle and crustal
%%% P velocities of 7.8 and 6.5 km/s (Ramachandran et al. 2006).
v_p = 7000;         % p-wave velocity [m/s]
%%% S-wave velocity calculated using v_p/v_s ratios from Audet et al. 2009. 
v_s = v_p/2.35;     % s-wave velocity [m/s] 
alpha = 0.2;        % Biot parameter [-] 
eta = 1e15;         % bulk viscosity [Pa*s]
mu = 2e-4;          % fluid viscosity [Pa*s]
nu = 0.3894;        % Poisson's ratio, average from Audet 2009; Peacock.
phi_i = 0.02;       % initial porosity [-] .04
phi_0 = .015;       % porosity cutoff [-] .015
rho_f = 1e3;        % fluid density [kg/m^3]
rho_s = 3e3;        % rock density [kg/m^3]

%%% Mechanical parameters. Computed quantaties.
rho_b = (1-phi_i)*rho_s + phi_i.*rho_f;     % bulk density [kg/m^3]
G = rho_b.*v_s^2;                           % shear modulus of rock [Pa]
K_u = 2*(1+nu).*G/(3*(1-2*nu));             % undrained bulk modulus [Pa]
%%% Drained bulk modulus [Pa]
K = (1./(2*phi_i)).*((K_f+K_u).*phi_i - alpha*(1+phi_i)...
    .*K_f + (4*(alpha-1).*(phi_i-alpha).*phi_i.*K_f.*K_u +...
    (phi_i.*(K_f+K_u)-alpha*(1+phi_i).*K_f).^2).^(1/2));
%%% Combined elastic parameter [Pa]
c_m = (K_u - K).*(K + 4*G/3)./(alpha^2*(K_u + 4*G/3));
%gamma = (K_u - K)./(alpha*(K_u + 4*G/3));   % Loading efficiency [-]

%%% Dehydration parameters for Lz + Brc = 2Fo + 3H2O (Wegner and Ernst,
%%% Ague et al.).
A = 201;            % surface area of rate limiting mineral [1/m] Connolly
c0 = 1.29e-20;      % reaction constant at the ref. temperature [mol/(m^2 s)
                    % (cal/mol)^n_r] Ague
E_a = 2e4;          % reaction activation energy [cal/mol] Ague
L = 500e3;          % latent heat of dehydration reaction [J/kg] Brantut
M_H2O = .018;       % molar mass of water [kg/mol]
n_r = 3.64;         % reaction order constant [-] Ague
R = 1.987;          % gas constant [cal/(K mol)]
T0 = 633;           % reference temperature [K] Ague
DoubleSnake = 0.1;  % plastic compaction parameter [-] Brantut

%%% Gibbs free energy parameters for reaction Atg + 20Brc = 34Fo + 51H2O
%%% (from supcrt).
a0 = 65945.7007021266;
a1 = -415.847466841169;
a2 = 0.654522687965514;
dT = 58;
a0 = a0 - a1*dT - 2*a2*dT + a2*dT^2;    %shifted value
a1 = a1-2*a2*dT;                        %shifted value
g0 = -643.932375104237;

%%% Intial conditions. The last 100 nodes do not contain dehydrating
%%% material.
N = 10000;      % Number of dehydrating nodes in the column
NE = 0.1*N;     % Number of non-dehydrating nodes at the base of the column
dz = 1000/N;    % distance between nodes [m]
N = 2*N + NE;   % total number of nodes in the column
Depth = (dz:dz:dz*N)';      % depth referenced to the top of the column [m]
    
%%% This function computes the initial conditions in the column.
    [pressure_i,p_lith,p_hydro,temperature_i,dTempdt,dsigmadt,dp_lithdz]...
        = Initial(phi_i,N,Depth,28.345e3,rho_s,rho_f);
%%% Initial porosity is uniform throughout the column.    
    phi_i = phi_i*ones(N,1);
%%% Initial released water content.
    m_di = zeros(N,1);
    
%%% Construct the initial conditions vector and solve.
if nargin < 2
    Vars_i = [pressure_i;temperature_i;phi_i;m_di];
    Time = 31536000*linspace(0,120000,2000)';           % Time vector [s]
else
    Vars_i = Vars_in;
    Time = 31536000*linspace(T_in,T_end,1401)';
end

%%% Differentiation matrices for finite differences.
%%% Set v = -1 to indicate that waves are moving towards shallower depths 
%%% (i.e. in negative z direction).
    v = -1;
%%% Three point upwind difference for first derivative    .
    D1 = three_point_upwind_uni_D1(Depth(1,1),Depth(N,1),N,v);
%%% Three point centered difference for second derivative.
    DD2 = three_point_centered_uni_D2(Depth(1,1),Depth(N,1),N);
    
%%% Generate the sparsity pattern for amazing speed up!
    S = SparsityPattern(N);
    options = odeset('Vectorized','on','JPattern',S,'RelTol',1e-10,'AbsTol',...
        [1e-12*ones(N,1); 1e-12*ones(N,1); 1e-12*ones(N,1); 1e-12*ones(N,1)]);
%%% Solve the system.    
    [Time,Vars_f] = ode15s(@PDE,Time,Vars_i,options);
    Time = Time/31536000;
    
%%% Deal with the output.    
    Pressure = Vars_f(:,1:N)';
%%% Pressure offset from initial condition.    
    PressureOffset = Pressure - repmat(pressure_i,1,size(Time,1));
    Load = dsigmadt*Time*31536000; % Loading from ongoing subduction.
%%% Pressure offset taking the loading into account.
    PressureOffset2 = PressureOffset - repmat(Load',numel(Depth),1);
%%% Temperature and porosity.    
    Temperature = Vars_f(:,N+1:2*N)';    
    Porosity = Vars_f(:,2*N+1:3*N)';
%%% Released water content.
    M_D = Vars_f(:,3*N+1:4*N)';
%%% Total fluid pressure.    
    P_tot = Pressure + repmat(p_hydro,1,size(Pressure,2));
%%% Effective pressure
    P_eff = repmat(p_lith,1,size(Time,1)) + repmat(Load',N,1)...
        - (Pressure + repmat(p_hydro,1,size(Time,1)));
%%% Gibbs free energy and reaction indexes.  
    DELTAG = -g0*((1/(2*a2))*(-a1 + (a1^2-4*a2*(a0-P_tot/1e5)).^(1/2))...
        - (Temperature-273));
    SIndex2 = DELTAG >= 0;  %reaction does not proceed if true.
%%% Kinetic dehydration equation.    
    Dehyd = (1/m_d0)*M_H2O*A*c0*(m_d0 - M_D).*(abs(DELTAG).^n_r)...
        .*exp((-E_a/R)*(1./Temperature - 1./T0));
    Dehyd(SIndex2) = 0;

function dVarsdt = PDE(time,Vars,~)
    p = Vars(1:N,:);
    T = Vars(N+1:2*N,:);
    phi = Vars(2*N+1:3*N,:);
    m_d = Vars(3*N+1:4*N,:);

%%% Calculate derivatives using differential matrices.
    dpdz = D1*p;
    dpdzz = DD2*p;
    dphidz = D1*phi;
    dTdzz = DD2*T;
    
%%% Dehydration equation.    
%%% Total fluid pressure.
    p_tot = repmat(p_hydro,1,size(p,2)) + p;
%%% Gibbs free energy of the dehydration reaction.    
    DeltaG = -g0*((1/(2*a2))*(-a1 + (a1^2 - 4*a2*(a0 - p_tot/1e5)).^(1/2))...
        - (T - 273));                                         %[cal/mol]
%%% Reaction does not proceed of DeltaG >= 0.    
    sIndex2 = DeltaG >= 0;   
%%% Calculate the reaction rate.    
    dm_fdt = (1/m_d0)*M_H2O*A*c0*(m_d0-m_d).*(abs(DeltaG).^n_r)...
        .*exp((-E_a/R)*(1./T - 1./T0));                     %[kg/(m^3 s)]
%%% Nodes with DeltaG >= 0 and bottom NE nodes are not dehydrating.    
    dm_fdt(N-NE+1:N,:) = zeros(NE,size(p,2));
    dm_fdt(sIndex2) = 0;
    
%%% Linear viscous relation for porosity change dPHI =
%%% -(3/4)*phi*p_eff/eta.
    p_eff = repmat(p_lith + dsigmadt*time,1,size(p,2)) - p_tot;
    dPHIdt = (DoubleSnake/rho_f)*dm_fdt - (3/4)*phi.*p_eff/eta;

%%% Pressure difference equation.
%%% Interior nodes.       
    dpdt(1:N-NE,:) = (k_i*c_m./mu).*(phi(1:N-NE,:)/phi_0).^n...
        .*((n./phi(1:N-NE,:)).*dphidz(1:N-NE,:).*dpdz(1:N-NE,:)...
       + dpdzz(1:N-NE,:)) + c_m.*(dm_fdt(1:N-NE,:)/rho_f...
       - dPHIdt(1:N-NE,:)) + dsigmadt; 
    dpdt(N-NE+1:N,:) = (k_i*c_m./mu).*(phi(N-NE+1:N,:)/phi_0).^n...
        .*((n./phi(N-NE+1:N,:)).*dphidz(N-NE+1:N,:).*dpdz(N-NE+1:N,:)...
       + dpdzz(N-NE+1:N,:)) + dsigmadt;
%%% Upper, absorbing boundary condition.
    dpdt(1,:) = (k_i*c_m./mu).*(phi(1,:)/phi_0).^n...
       .*(n./phi(1,:)).*dphidz(1,:).*dpdz(1,:)...
       + (c_m/eta)*(3/4)*phi(1,:).*p_eff(1,:) + dsigmadt;
%%% Bottom boundary. Lithostatic gradient constant pressure condition.
    pN = pressure_i(N,1) + dz*dp_lithdz + dsigmadt*time;
    dpdzzN = (p(N-1,:) - 2*p(N,:) + pN)/dz^2;  
    dpdt(N,:) = (k_i*c_m./mu).*(phi(N,:)/phi_0).^n.*((n./phi(N,:)).*...
        dphidz(N,:).*dpdz(N,:) + dpdzzN) + c_m.*(dm_fdt(N,:)/rho_f...
        - dPHIdt(N,:)) + dsigmadt; 

%%% Porosity equation.    
    dphidt(1:N-NE,:) = -(alpha./(K+4*G/3)).*dsigmadt + (1./c_m...
        - phi(1:N-NE,:)./K_f).*dpdt(1:N-NE,:) + dPHIdt(1:N-NE,:);
    dphidt(N-NE+1:N,:) = -(alpha./(K+4*G/3)).*dsigmadt + (1./c_m...
        - phi(N-NE+1:N,:)./K_f).*dpdt(N-NE+1:N,:);

%%% Temperature equation.
%%% Only heat conduction.
    dTdt = (1/C_s)*(K_T*dTdzz - L*dm_fdt);
%%% Upper boundary.
    dTdt(1,:) = -(L/C_s)*dm_fdt(1,:);    
%%% Constant temperature increase at bottom boundary.
    dTdt(N,:) = dTempdt;
%%% Assemble all of the variables.    
    dVarsdt = [dpdt;dTdt;dphidt;dm_fdt];
end %PDE
end

function [Pressure_i,Pressure_0,Hydrostatic,Temp,dTempdt,dsigmadt,dp_lithdz]...
    = Initial(phi_i,N,Depth,Depth_0,rho_s,rho_f)
%%% This function computes the initial conditions in the column based on
%%% the starting depth of the simulation and the total depth. And also
%%% computes some parameter values related to ongoing subduction.

%%% INPUT
% phi_i                         Initial porosity in the column [-]
% N                             Number of nodes in the column
% Depth                         Depth referenced to the top of the column [m]
% Depth_0                       Depth of the top of the column referenced
%                               to the surface [m]
% rho_f                         fluid density [kg/m^3]
% rho_s                         rock density [kg/m^3]

%%% OUTPUT
% Pressure_i                    Excess fluid pressure in column [Pa]
% Pressure_0                    Total fluid pressure in column [Pa]
% Hydrostatic                   Hydrostatic fluid pressure in column [Pa]
% Temp                          Temperature in column [K]
% dTempdt                       Rate of temperature increase at the base of
%                               the column due to subduction [K/s]
% dsigmadt                      Rate of increase in the overburden due to
%                               subduction [Pa/s]
% d_plithdz                     Lithostatic overpressure gradient [Pa/m]


%%% Some examples:
%%% Initial pressure and temperature profiles. 
%%% For Depth = 1000 m:
%%% 0 nodes dehydrating: 32.19355e3
%%% 1 node dehydrating: 32.1962e3
%%% About 12 m dehydrating: 32.24e3
%%% About 50 nodes dehydrating: 32.4e3

%%% For Depth  = 2000 m
%%% About 50 m dehydrating: 28.5e3
%%% About 12 m dehydrating: 28.345e3, gives same DeltaG_i profile as
%%% 32.24e3 for 1000 m column.

%%% For Depth = 3000
%%% About 12 m dehydrating: 24.383e3, gives same DeltaG_i profile as
%%% 32.24e3 for 1000 m column.

%%% Plate interface for distances far from the trench. The coefficients of 
%%% the plate interface polynomial are found from a linear fit to the plate
%%% interface depth for depths ~20-30 km
    PlateInt = [0.2709 -9.8967e3]; %(meters)%
    Distance_0 = (Depth_0 - PlateInt(2))/PlateInt(1);
    
%%% Geometrical Parameters for Temperature and LoadingConstant calculation.
    betaT = atan(PlateInt(1));
    v_plate = 0.037/(365*24*3600); % plate convergence velocity [m/s]
    v_x = v_plate*cos(betaT);
    StressFit = [8693.11664571161 -365056544.708814];
    dsigmadt = StressFit(1)*cos(betaT)*v_plate;
   
%%% Temperature field estimated from Peacock 2009. TempDepth is the linear
%%% increase with subduction depth.
    TempDepth = [200/30000 (400*50000-600*20000)/30000];
    dTdz = 500/25000;
    dTempdt = TempDepth(1)*PlateInt(1)*v_x;
    Temp = TempDepth(1).*(PlateInt(1).*Distance_0+PlateInt(2))+TempDepth(2)...
        + dTdz.*Depth + 273;
%%% Calculate the stress fit and initial overburden. The stress coefficients
%%% are found from a linear fit the stress estimated at the plate interface
%%% from the Trehu velocity model, for the same range as the plate interface
%%% polynomial.
    StressFit = [8693.11664571161 -365056544.708814]; %(Pa%
    Overburden = polyval(StressFit,Distance_0);      
%%% Set the initial porosity.        
    Porosity = phi_i.*ones(N,1);
    
%%% Calculate the lithostatic stress in the column due to the self-weight,
%%% assuming an initial exess pore pressure, Pressure_0, in the column.
    Lithostatic = zeros(N,1);
    BulkDensity = rho_s - (rho_s-rho_f).*Porosity; 
    Lithostatic(1,1) = 9.81*(Depth(1))*BulkDensity(1);
        for i = 2:N
            Lithostatic(i,1) = 9.81*(Depth(i) - Depth(i-1))*BulkDensity(i)...
                + Lithostatic(i-1);                
        end    
%%% Hydrostatic pressure with reference datum set to the top of the column.    
    Hydrostatic = 9.81*rho_f*(Depth+Depth_0);
%%% This is the total pressure, the term proportional to the overburden is
%%% the initial excess pressure. For lithostatic pore pressure, the EXCESS
%%% pore pressure P_e = Overburden + Lithostatic - Hydrostatic; thus the 
%%% TOTAL pore pressure is P_t = Overburden + Lithostatic, and the
%%% effective stress will then be zero everywhere.
%%% Pressure_0 is the TOTAL initial pore pressure.
    Pressure_0 = (Overburden).*ones(N,1) + Lithostatic;
    Pressure_i = Pressure_0 - Hydrostatic;
%%% Calculate the lithostatic overpressure gradient.    
    dp_lithdz = 9.81*(BulkDensity - rho_f);
    dp_lithdz = mean(dp_lithdz);
end %Initial

function S = SparsityPattern(N)
%%% Numerical set up. Create the JPattern sparsity matrix based on the
%%% finite-difference equations.
    iD = (1:N);
    iSup = (1:N-1);
    iSupSup = (1:N-2);
%%% Pressure section.
    P_p = sparse(iD,iD,ones(N,1)) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSup+1,iSup,ones(N-1,1),N,N)...
        + sparse(iSupSup,iSupSup+2,ones(N-2,1),N,N);
    P_p(N,N-2) = 1;
    P_p(N,N-3) = 1;
    P_T = sparse(iD,iD,ones(N,1));
    P_phi = sparse(iD,iD,ones(N,1),N,N) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSupSup,iSupSup+2,ones(N-2,1),N,N);
    P_phi(N-1,N-2) = 1;
    P_phi(N,N-1) = 1;
    P_phi(N,N-2) = 1;
    P_m = P_T;
    Jac_P = [P_p P_T P_phi P_m];
%%% Temperature equation. 
    T_p = sparse(iD,iD,[ones(N-1,1);0]) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSup+1,iSup,[ones(N-2,1);0],N,N);
    T_T = T_p;
    T_phi = P_phi;
    T_m = sparse(iD,iD,[ones(N-1,1);0]);
    Jac_T = [T_p T_T T_phi T_m];
%%% Porosity section is equal to the pressure section.
%%% Dehydration section.
    Jac_m = [P_T P_T sparse(N,N) P_T];
%%% Combine all the sections to form the sparsity pattern.
    S = [Jac_P; Jac_T; Jac_P; Jac_m];
end %SparsityPattern