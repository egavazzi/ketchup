%GENERAL PARAMETERS
dt = 2.5e-5;                % timestep
Niter=16000;                % number of iterations to perform
dump_period_fields = 40;    % number of iterations between field dumps
fields_per_file = 5;        % number of field dumps to gather in each file
dump_period_distr = 800;    % number of iterations between distribution dumps
dump_period_distr_1v = 400; % number of iterations between f(z,vz) dumps
dump_start = 1;          % iteration at which dumping shall commence
shift_test_period = 10;  % number of time steps between testing for the
                         % need to shift and between updates of the voltage 
                         % across the series resistance
resistance = 0.0e9;      % series resistance in ohms for one square metre
                         % of fluxtube cross section at the left hand side
Nz=875;                  % number of grid points in space
zmin=0.0;                % left hand boundary
zmax=5.5e7;              % right hand boundary
Nspecies=4;              % number of species
const_a = 400.0e0;       % epsilon_r=max((omegap*dt*const_a)**2,1.0d0)
BC_Poisson=1;            % Boundary condition for Poisson's equation.
                         % 1: V(right hand side) = voltage
                         % 2: E(left hand side) = E0
                         % 3: E(right hand side) = E0
voltage = 3.0e3;         % voltage across the system (potential of right
                         % hand side)
initialiser = 1;         % your choice of initialisation routine
                         % 1: original version
E0=0.0e0;                % Field at left hand boundary when BC_Poisson=2.
startfromdumpfile='yes'; % if 'yes' the distribution is read from the
                         % dumps directory
dump_period_dump=800;    % number of iterations between dumps you can
                         % start from
exit_after_dump='yes';   % whether execution shall stop after the first dump
transffilename='polestestcase1.dat'; % Name of transform file
voltagefilename='no file';    % Name of voltage control file. If
                              % missing or file non-existent a
                              % constant voltage=voltage is used.
%END

%SPECIES 1: magnetospheric electrons
Nvz=114;               % number of grid points in velocity
vzmin =   -6.84e7;      % minimum velocity - vz0
vzmax =   6.84e7;       % maximum velocity - vz0
vzshifting = 'off';
Nmu = 50;              % number of grid points in magnetic moment
mumin = 0;             % minimum magnetic moment
mumax = 7.0e-9;        % maximum magnetic moment
muexp = 4; % exponent for nonuniform mu-grid mu~([1:Nmu]/Nmu)^muexp
mass = 9.10938215e-31; % if mass is NaN, an infinite mass neutralising
                       % background is assumed
charge = -1.602176487e-19;
relativistic = 'yes';  % if 'yes', the particle mass is computed
                       % according to the theory of relativity.
n0 = 3.0e5;            % initial density maximum
vfromfile = 'no';      % If 'yes', load initial vz from file. The global vz0
                       % will still be the reference for velocity shifting.
vz0 = 0.0e0;           % drift velocity
Tfromfile = 'no';      % If 'yes', load T from file and ignore kTz and kTp.
kTz = 5.0e2;    % voltage equivalent to temperature kTz/e. f~exp(mv^2/(2kTz))
kTp = 5.0e2;    % voltage equivalent to temperature kTp/e. f~exp(mv^2/(2kTp))
n0L = 3.0e5;    % density at left hand boundary
vz0L = 0.0e6;   % drift velocity at left hand boundary
kTzL = 5.0e2;   % parallel temperature at left hand boundary
kTpL = 5.0e2;   % perpendicular temperature at left hand boundary
lossconeL = 'no'; % if 'yes' the left hand boundary loss cone is empty.
n0R = 0.0e9;      % density at right hand boundary
vz0R = 0.0e5;     % drift velocity at right hand boundary
kTzR = 5.0e2;     % parallel temperature at right hand boundary
kTpR = 5.0e2;     % perpendicular temperature at right hand boundary
lossconeR = 'no'; % if 'yes' the right hand boundary loss cone is empty.
%END

%SPECIES 2: magnetospheric protons
Nvz=100;               % number of grid points in velocity
vzmin =   -6.0e6;      % minimum velocity - vz0
vzmax =    6.0e6;      % maximum velocity - vz0
vzshifting = 'off';
Nmu = 50;              % number of grid points in magnetic moment
mumin = 0;             % minimum magnetic moment
mumax = 3.5e-8;        % maximum magnetic moment
muexp = 4; % exponent for nonuniform mu-grid mu~([1:Nmu]/Nmu)^muexp
mass = 1.672621637e-27;% if mass is NaN, an infinite mass neutralising
                       % background is assumed
%mass = NaN;
charge = 1.602176487e-19;
relativistic = 'no';   % if 'yes', the particle mass is computed
                       % according to the theory of relativity.
n0 = 3.0e5;            % initial density maximum
vz0 = 0.0e0;           % drift velocity
Tfromfile = 'no';      % if 'yes', load T from file and ignore kTz and kTp.
kTz = 2.5e3;    % voltage equivalent to temperature kTz/e. f~exp(mv^2/(2kT))
kTp = 2.5e3;    % voltage equivalent to temperature kTp/e. f~exp(mv^2/(2kT))
n0L = 3.0e5;    % density at left hand boundary
vz0L = 0.0;     % drift velocity at left hand boundary
kTzL = 2.5e3;   % parallel temperature at left hand boundary
kTpL = 2.5e3;   % perpendicular temperature at left hand boundary
lossconeL = 'no'; % if 'yes' the left hand boundary loss cone is empty.
n0R = 0.0e9;      % density at right hand boundary
vz0R = 0.0e0;     % drift velocity at right hand boundary
kTzR = 2.5e3;     % parallel temperature at right hand boundary
kTpR = 2.5e3;     % perpendicular temperature at right hand boundary
lossconeR = 'no'; % if 'yes' the right hand boundary loss cone is empty.
%END

%SPECIES 3: ionospheric electrons
Nvz=100;               % number of grid points in velocity
vzmin =   -1.0e7;      % minimum velocity - vz0
vzmax =   1.0e7;       % maximum velocity - vz0
vzshifting = 'off';
Nmu = 20;              % number of grid points in magnetic moment
mumin = 0;             % minimum magnetic moment
mumax = 2.5e-14;       % maximum magnetic moment
muexp = 4;             % exponent for nonuniform mu-grid mu~([1:Nmu]/Nmu)^muexp
mass = 9.10938215e-31; % if mass is NaN, an infinite mass neutralising
                       % background is assumed
charge = -1.602176487e-19;
relativistic = 'no';   % if 'yes', the particle mass is computed
                       % according to the theory of relativity.
n0 = 1.542857e6;       % initial density maximum
vfromfile = 'no';      % If 'yes', load initial vz from file. The global vz0
                       % will still be the reference for velocity shifting.
vz0 = 0.0e0;           % drift velocity
Tfromfile = 'no';      % If 'yes', load T from file and ignore kTz and kTp.
kTz = 1.0e0;    % voltage equivalent to temperature kTz/e. f~exp(mv^2/(2kTz))
kTp = 1.0e0;    % voltage equivalent to temperature kTp/e. f~exp(mv^2/(2kTp))
n0L = 0.0e5;    % density at left hand boundary
vz0L = 0.0e6;   % drift velocity at left hand boundary
kTzL = 5.0e2;   % parallel temperature at left hand boundary
kTpL = 5.0e2;   % perpendicular temperature at left hand boundary
lossconeL = 'no'; % if 'yes' the left hand boundary loss cone is empty.
%n0R = 1.0e9;     % density at right hand boundary
n0R = 1.542857e6; % density at right hand boundary
vz0R = 0.0e5;     % drift velocity at right hand boundary
kTzR = 1.0e0;     % parallel temperature at right hand boundary
kTpR = 1.0e0;     % perpendicular temperature at right hand boundary
lossconeR = 'no'; % if 'yes' the right hand boundary loss cone is empty.
%END

%SPECIES 4: ionospheric protons
Nvz=100;               % number of grid points in velocity
vzmin =   -2.5e5;      % minimum velocity - vz0
vzmax =    2.5e5;      % maximum velocity - vz0
vzshifting = 'on';
Nmu = 20;              % number of grid points in magnetic moment
mumin = 0;             % minimum magnetic moment
mumax = 2.5e-14;       % maximum magnetic moment
muexp = 4;             % exponent for nonuniform mu-grid mu~([1:Nmu]/Nmu)^muexp
mass = 1.672621637e-27;% if mass is NaN, an infinite mass neutralising
                       % background is assumed
%mass = NaN;
charge = 1.602176487e-19;
relativistic = 'no';   % if 'yes', the particle mass is computed
                       % according to the theory of relativity.
%n0 = 1.0e9;           % initial density maximum
n0 = 1.542857e6;       % initial density maximum
vz0 = 0.0e0;           % drift velocity
Tfromfile = 'no';      % if 'yes', load T from file and ignore kTz and kTp.
kTz = 1.0e0;    % voltage equivalent to temperature kTz/e. f~exp(mv^2/(2kT))
kTp = 1.0e0;    % voltage equivalent to temperature kTp/e. f~exp(mv^2/(2kT))
n0L = 0.0e5;    % density at left hand boundary
vz0L = 0.0;     % drift velocity at left hand boundary
kTzL = 5.0e2;   % parallel temperature at left hand boundary
kTpL = 5.0e2;   % perpendicular temperature at left hand boundary
lossconeL = 'no';      % if 'yes' the left hand boundary loss cone is empty.
%n0R = 1.0e9;          % density at right hand boundary
n0R = 1.542857e6;      % density at right hand boundary
vz0R = 0.0e0;          % drift velocity at right hand boundary
kTzR = 1.0e0;          % parallel temperature at right hand boundary
kTpR = 1.0e0;          % perpendicular temperature at right hand boundary
lossconeR = 'no';      % if 'yes' the right hand boundary loss cone is empty.
%END
