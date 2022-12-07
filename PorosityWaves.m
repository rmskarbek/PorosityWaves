function [Time,Pressure,PressureOffset,Porosity,Depth,P_tot,P_eff]...
    = PorosityWaves
%%% This code solves the dimensionaless equations for pore pressure and
%%% porosity, for an initial porosity step, without reactions. There is an
%%% absorbing condition at the upper boundary. Pressure is constant at the
%%% lower boundary.

%%% Mechanical parameters. Constants.
g = 9.81;           % gravity [m/s^2]
n = 3;              % perm-porosity power law exponent [-]
phi_i = 0.02;       % initial porosity [-]
phi_1 = 0.8;        % magnitude of dimensionless porosity step.
rho_f = 1e3;        % fluid density [kg/m^3]
rho_s = 3e3;        % rock density [kg/m^3]

%%% Set your own dimensionless parameters.
    delta = 30;                             % viscous compaction length
    p_star = g*(rho_s - rho_f)*delta;       % pressure scale
    De_P = 1e-4;                            % pressure Deborah number
    De_PHI = 1e-4;                          % porosity Deborah number

%%% Intial conditions.
%%% Grid spacing must be some amount smaller than the compaction length,
%%% according to Spiegelman 1993b.
    dz = delta/(2*4*phi_1^(-3/2));
    N = round(200*delta/dz);            % total grid points
    Depth = (dz:dz:dz*N)';              % dimensional depth
%%% Initial pressure profiles.
    [pressure_i,p_lith,p_hydro,~,~,~,~]...
        = Initial(phi_i,N,Depth,32.1962e3,rho_s,rho_f);
%%% Dimensionless depth coordinates.
    dz = dz/delta;
    Depth = (dz:dz:dz*N)';
    
%%% Construct the initial conditions vector and solve.
    pressure_i = pressure_i/p_star;
    phi_I = ones(N,1);
%%% Initial porosity step.    
    N1 = round(4*N/5);
    phi_I(1:N1) = phi_1 + (1 - phi_1)*sech((Depth(1:N1) - Depth(N1))/2.5);
    Vars_i = [pressure_i;phi_I];
    Time = linspace(0,150,1000)';            % dimensionless time vector.

%%% Differentiation matrices for finite differences
    v = -1;
    D1 = three_point_upwind_uni_D1(Depth(1,1),Depth(N,1),N,v);
%%% Three point centered difference for second derivative
    DD2 = three_point_centered_uni_D2(Depth(1,1),Depth(N,1),N);
    
%%% Generate the sparsity pattern for amazing speed up!
    S = SparsityPattern(N);
    options = odeset('Vectorized','on','JPattern',S,'RelTol',1e-6,'AbsTol',...
        1e-9,'Stats','on');    
    [Time,Vars_f] = ode15s(@PDE,Time,Vars_i,options);
    
%%% Deal with the output.    
    Pressure = Vars_f(:,1:N)';
    PressureOffset = Pressure - repmat(pressure_i,1,size(Time,1));
    Porosity = Vars_f(:,N+1:2*N)';
%%% Total fluid pressure.    
    P_tot = Pressure + repmat(p_hydro/p_star,1,size(Pressure,2));
%%% Porosity change as a function of effective stress and mass transefer
%%% due to dehydration.
    P_eff = repmat(p_lith/p_star,1,size(Time,1))...
        - (Pressure + repmat(p_hydro/p_star,1,size(Time,1)));

function dVarsdt = PDE(~,Vars,~)
%%% These variables are dimensionless.    
    p = Vars(1:N,:);
    phi = Vars(N+1:2*N,:);
%%% Dimensionlesss compaction calculations.
    p_tot = repmat(p_hydro/p_star,1,size(p,2)) + p;
%%% Linear viscous relation for porosity change.
    p_eff = repmat(p_lith/p_star,1,size(p,2)) - p_tot;  
  
%%% First spatial derivative of pressure. Two point upstream for a 
%%% non-uniform grid. Set v = -1 to indicate that waves are moving towards
%%% shallower depths (i.e. in negative z direction).
    dpdz = D1*p;
    dpdzz = DD2*p;
    dphidz = D1*phi;

%%% Compute governing equations.
%%% Pressure equation and BCs.
    dpdzz(1,:) = (p(1,:) - 2*p(2,:) + p(3,:))/dz^2;
    dpdt = (1/De_P)*(n*phi.^(n-1).*dphidz.*dpdz + phi.^n.*dpdzz + phi.*p_eff);
%%% Upper, absorbing boundary condition.  
    dpdt(1,:) = (1/De_P)*(n*phi(1,:).^(n-1).*dphidz(1,:).*dpdz(1,:)...
        + phi(1,:).*p_eff(1,:));
%%% Lower boundary.    
    dpdt(N,:) = zeros(1,size(p,2));
   
%%% Porosity equation.    
    dphidt = (De_P - De_PHI*phi).*dpdt - phi.*p_eff;  
%%% Assemble all of the variables.       
    dVarsdt = [dpdt;dphidt];
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
    iD = (1:N);
    iSup = (1:N-1);
    iSupSup = (1:N-2);
%%% Numerical set up. Create the JPattern sparsity matrix.
%%% Pressure section.
%%% Dirichlet boundary at top and bottom.
    P_p = sparse(iD,iD,[ones(N-1,1);0],N,N) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSup+1,iSup,[ones(N-2,1);0],N,N)...
        + sparse(iSupSup,iSupSup+2,ones(N-2,1),N,N);
    P_p(1,4) = 1;
    P_phi = sparse(iD,iD,[ones(N-1,1);0],N,N) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSupSup,iSupSup+2,ones(N-2,1),N,N);
    P_phi(N-1,N-2) = 1;
    Jac_P = [P_p P_phi];
%%% Porosity section.
%%% Dirichlet boundaries.
    Phi_p = sparse(iD,iD,ones(N,1),N,N) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSup+1,iSup,[ones(N-2,1);0],N,N)...
        + sparse(iSupSup,iSupSup+2,ones(N-2,1),N,N);
    Phi_p(1,4) = 1;
    Phi_phi = sparse(iD,iD,[ones(N-1,1);0],N,N) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSupSup,iSupSup+2,ones(N-2,1),N,N);
    Phi_phi(N-1,N-2) = 1;
    Jac_phi = [Phi_p Phi_phi];
    
%%% Full J pattern.
    S = [Jac_P; Jac_phi];

end %SparsityPattern