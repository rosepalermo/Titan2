% calculate change in area for circle to make sure our model has the proper
% units and for a grid bias test.

% accompanies slideshow "Analytical error of model 0621"

% for the paper, account for "erosion" or strength loss at the shoreline
% that would be partial erosion of the cell. --- NAH --- 

% notes for Rose - dividing by 1+sqrt(2) in addidshoreline -- this is
% correct
% when don't divide by sqrt(2), change in area is about double analytical
% (1.926x)-->darea_analytical = 1.2; darea_model = 2.4
% when divide by sqrt(2), darea_uniform_analytical = 1.2. model = 1.4

% now notes on passing damage - minimizing the error with mean or max
% when max, danalytical = 13.06 and model = 15.66. when mean, model =
% 15.345. So mean is slightly closer

%% UNIFORM -- corners (con8 = 0)

load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/circle_analytical_error/uniform_con8_0_tf1000_d1sqrt2_mean.mat')

% ANALYTICAL
% initial area of cirlce is 25pi
r0 = 5;
Area0 = pi*(r0^2);
k_uniform = p.Kuniform*(1/(p.So./p.dxo));
t = 4000;
Areat_uniform = pi*(r0+k_uniform*t).^2;
t_erode_1cell = (p.dx)/(k_uniform);
dArea_uniform_analytical = Areat_uniform(end) - Area0;

% MODELED
uniform_init = ~(output(:,:,1));
uniform_final = g.U<1;
uniform_init_area = length(find(uniform_init))./((length(uniform_init)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
uniform_final_area = length(find(uniform_final))./((length(uniform_final)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
% sl_lessthanstrength = g.Strength(g.Strength<p.strength);
% sl_lessthanstrengthandgreaterthan0 = sl_lessthanstrength(sl_lessthanstrength>0);
% strengthloss_shoreline = sum(p.strength-sl_lessthanstrengthandgreaterthan0);
% equivalent_cells = strengthloss_shoreline/p.strength;
% equivalent_cells_area = equivalent_cells./((length(uniform_final)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
dArea_uniform_model = uniform_final_area - uniform_init_area;


%% WAVE -- corners (con8 = 0)
% This is more complicated than this.. Assuming that total damage per time
% step is Kwave*dt.. But there is more in there. So 

load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/circle_analytical_error/wave_con8_0_tf1000.mat')

% ANALYTICAL
r0 = 5;
Area0 = pi*(r0^2);
k_wave = p.Kwave./Area0*(1/(p.So./p.dxo));
t = 1000;
Areat_wave = pi*((r0/(1-r0*k_wave*3/4*pi*t)).^2);
dArea_wave_analytical = Areat_wave(end) - Area0;

% MODELED
wave_init = ~(output(:,:,1));
wave_final = g.U<1;
wave_init_area = length(find(wave_init))./((length(wave_init)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
wave_final_area = length(find(wave_final))./((length(wave_final)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
dArea_wave_model = wave_final_area - wave_init_area;

%% UNIFORM -- no corners (con8 = 1)
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/circle_analytical_error/uniform_con8_1_tf1000_dsqrt2.mat')

% ANALYTICAL
% initial area of cirlce is 25pi
r0 = 5;
Area0 = pi*(r0^2);
k_uniform = p.Kuniform*(1/(p.So./p.dxo));
t = 1000;
Areat_uniform = pi*(r0+k_uniform*t).^2;
t_erode_1cell = (r0+p.dx-r0)/(k_uniform);
dArea_uniform_analytical = Areat_uniform(end) - Area0;

% MODELED
uniform_init = ~(output(:,:,1));
uniform_final = g.U<1;
uniform_init_area = length(find(uniform_init))./((length(uniform_init)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
uniform_final_area = length(find(uniform_final))./((length(uniform_final)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
dArea_uniform_model = uniform_final_area - uniform_init_area;

%% WAVE -- no corners (con8 = 1)
% This is more complicated than this.. Assuming that total damage per time
% step is Kwave*dt.. But there is more in there. So 

load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/circle_analytical_error/wave_con8_1_tf1000.mat')

% ANALYTICAL
r0 = 5;
Area0 = pi*(r0^2);
k_wave = p.Kwave./Area0*(1/(p.So./p.dxo));
t = 1000;
Areat_wave = pi*((r0/(1-r0*k_wave*3/4*pi*t)).^2);
dArea_wave_analytical = Areat_wave(end) - Area0;

% MODELED
wave_init = ~(output(:,:,1));
wave_final = g.U<1;
wave_init_area = length(find(wave_init))./((length(wave_init)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
wave_final_area = length(find(wave_final))./((length(wave_final)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
dArea_wave_model = wave_final_area - wave_init_area;