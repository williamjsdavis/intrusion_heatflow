%% Intrusion heat flow
% William Davis 30/06/16
% Pascal and Olesen 2009
clc, clearvars, close all

%% Variables
x_length = 500E3; % m
y_length = 140E3; % m
crust_thickness = 40E3; % m
diapir_length = 100E3; % m
scale = 0.0001; % Points per m
nt = 200; % Number of timesteps
eps = 0.5; % Stability limit

%% Paper variables
% Crustal heat generation
A0 = 0.35; % milli watts per meter cubed

% Temperature of the diapir
T_d = 1500; % Degrees C

% Thermal diffusivity
tDifCrust = 0.9E-6; % Crust % Thermal diffusivity, m^2 s^-1
tDifMant = 1E-6; % Mantle % Thermal diffusivity, m^2 s^-1

% Thermal conductivity
tConCrust = 2.5; % Crust
tConMant = 3.5; % Mantle

%%
nptsx = x_length*scale; % Calculates grid spacing for rectangular boxes
nptsy = y_length*scale;
nx = nptsx + 1; % Number of x and y nodes
ny = nptsy + 1;
nnodes = nx*ny; % Total number of nodes
dx = x_length/(nx-1);  % Calculate x space increment in km
dy = y_length/(ny-1);  % Calculate y space increment in km
kmSca = 1000; % km to m scale


%% Setting up grid
thCrust = ceil(ny*crust_thickness/y_length); % Thickness of crust
thMant = ny - thCrust; % Thickness of mantle
lenDia = ceil(nx*diapir_length/x_length); % Length of diapir
lenMant = (nx - lenDia)/2;

%% Assigning values

% Thermal diffusivity and conductivity grids
K = ones(ny,nx); % Diffusivity
k = ones(ny,nx); % Conductivity

K(1:thCrust,:) = tDifCrust;
K(thCrust+1:end,:) = tDifMant;
k(1:thCrust,:) = tConCrust;
k(thCrust+1:end,:) = tConMant;

K_max = max(K(:)); % Maximum diffusivity
d_t = eps*dx*dy/K_max; % Calculate time step from stability criterion 
                                   %(dimensional time)

dt = d_t/(365.25*24*60*60); % In years
dtMa = dt/1000000; % In Ma

% Setting up grid
[X,Y] = meshgrid((0:dx:x_length),(0:dy:y_length));

% Set up initial T
T_gradient = 0:T_d/(ny-1):T_d;
T0 = repmat(T_gradient',1,nx); % Set up crust and mantle gradient
T0(thCrust+1:end,lenMant:lenMant+lenDia) = T_d;
T = T0; T_new = T0; T_old = T0; % Initialization of temperature fields
T_mat = zeros(ny,nx,nt+1); T_mat(:,:,1) = T0; % Saved T matrix

%% QSM vectorised method

% Lookup number matrix
number_mat = zeros(size(X));
number_mat(:) = 1:numel(X);

% Diffusivity coefficient (s m^-2)
sx = d_t/(dx*dx);
sy = d_t/(dy*dy);

% Inbetween nodes for K
K_center = K(2:ny-1,2:nx-1);
K_north = K(1:ny-2,2:nx-1);
K_south = K(3:ny,2:nx-1);
K_east = K(2:ny-1,3:nx);
K_west = K(2:ny-1,1:nx-2);
Kxe_mat = (K_east+K_center)/2; % East
Kxw_mat = (K_west+K_center)/2; % West
Kys_mat = (K_south+K_center)/2; % South
Kyn_mat = (K_north+K_center)/2; % North

center_mat = number_mat(2:ny-1,2:nx-1);
north_mat = number_mat(1:ny-2,2:nx-1);
south_mat = number_mat(3:ny,2:nx-1);
east_mat = number_mat(2:ny-1,3:nx);
west_mat = number_mat(2:ny-1,1:nx-2);

% Edges
north_nodes = ny+1:ny:(nnodes-ny);
south_nodes = 2*ny:ny:(nnodes-ny);
east_nodes = nnodes-ny+2:nnodes-1;
west_nodes = 2:ny-1;
% Corners
nw_node = 1;
ne_node = nnodes-ny+1;
sw_node = ny;
se_node = nnodes;

ivec = ones(5*nnodes,1)*NaN;
jvec = ones(5*nnodes,1)*NaN;
vvec = ones(5*nnodes,1)*NaN;

ivec(5*center_mat(:)-4) = center_mat(:);
ivec(5*center_mat(:)-3) = center_mat(:);
ivec(5*center_mat(:)-2) = center_mat(:);
ivec(5*center_mat(:)-1) = center_mat(:);
ivec(5*center_mat(:)) = center_mat(:);
jvec(5*center_mat(:)-4) = south_mat(:);
jvec(5*center_mat(:)-3) = north_mat(:);
jvec(5*center_mat(:)-2) = center_mat(:);
jvec(5*center_mat(:)-1) = east_mat(:);
jvec(5*center_mat(:)) = west_mat(:);
vvec(5*center_mat(:)-4) = -sy*Kys_mat;
vvec(5*center_mat(:)-3) = -sy*Kyn_mat;
vvec(5*center_mat(:)-2) = (Kxe_mat+Kxw_mat)*sx+(Kys_mat+Kyn_mat)*sy+1;
vvec(5*center_mat(:)-1) = -sx*Kxe_mat;
vvec(5*center_mat(:)) = -sx*Kxw_mat;

% Corners
ivec(nw_node*5) = nw_node; % NW
jvec(nw_node*5) = nw_node;
vvec(nw_node*5) = 1;
ivec(ne_node*5) = ne_node; % NE
jvec(ne_node*5) = ne_node;
vvec(ne_node*5) = 1;
ivec(sw_node*5) = sw_node; % SW
jvec(sw_node*5) = sw_node;
vvec(sw_node*5) = 1;
ivec(se_node*5) = se_node; % SE
jvec(se_node*5) = se_node;
vvec(se_node*5) = 1;

% Edges
ivec(west_nodes*5) = west_nodes; % West boundary
jvec(west_nodes*5) = west_nodes;
vvec(west_nodes*5) = 1;
ivec(east_nodes*5) = east_nodes; % East boundary
jvec(east_nodes*5) = east_nodes;
vvec(east_nodes*5) = 1;

ivec(north_nodes*5) = north_nodes; % North boundary
jvec(north_nodes*5) = north_nodes;
vvec(north_nodes*5) = 1;
ivec(south_nodes*5) = south_nodes; % South boundary
jvec(south_nodes*5) = south_nodes;
vvec(south_nodes*5) = 1;

% Final calcs
ivec(isnan(ivec)) = [];
jvec(isnan(jvec)) = [];
vvec(isnan(vvec)) = [];
qsm = sparse(ivec,jvec,vvec);

%% Time stepping
for i = 1:nt 
   T_new = qsm\T(:);
   T_old = T;
   T(:) = T_new(:); % Change times around
   T_mat(:,:,i+1) = T; % Saving T matrix
   
   figure(1)
   contourf(X/1000,Y/1000,T,'LineStyle','none','LevelStep',100) % 1600/9
   %c = colorbar('southoutside','Limits',[0,1600]);
   c.Label.String = 'Temp. (\circC)';
   colormap(jet)
   ylabel('Depth, km')
   axis ij
   axis equal
   
   line([0,x_length/kmSca],[crust_thickness/kmSca,crust_thickness/kmSca],'Color','k',...
       'LineWidth',1)
   line([(x_length-diapir_length)/(2*kmSca),(x_length-diapir_length)/(2*kmSca)],...
       [crust_thickness/kmSca,y_length/kmSca],'Color','k',...
       'LineWidth',1)
   line([(x_length+diapir_length)/(2*kmSca),(x_length+diapir_length)/(2*kmSca)],...
       [crust_thickness/kmSca,y_length/kmSca],'Color','k',...
       'LineWidth',1)
   
   timestep = num2str(i,'%03.0f');
   ageMa = num2str(round(i*dtMa,1),'%04.1f');
   title([' FD, time step: ',timestep,'     Time(Ma): ',ageMa])
   xlabel ('Distance, km'), ylabel('Depth, km')
   drawnow
   
end
save('T_mat')