% calculate change in area for circle grid bias test for uniform and wave
% erosion

% accompanies slideshow "Analytical error of model"

folder = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/grid bias testing/circletest/';

%% UNIFORM
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/grid bias testing/circletest/uniform_Kc0_0001.mat')

% ANALYTICAL
% initial area of cirlce is 25pi
Area0 = 25*pi;
r0 = 5;
k_uniform = p.Kcoast;
t = 1;
Areat_uniform = pi*(r0+k_uniform*p.dx*t).^2;
dArea_uniform_analytical = Areat_uniform - Area0;

% MODELED
uniform_init = ~(output(:,:,1));
uniform_final = g.U<1;
uniform_init_area = length(find(uniform_init))./((length(uniform_init)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
uniform_final_area = length(find(uniform_final))./((length(uniform_final)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
dArea_uniform_model = uniform_final_area - uniform_init_area;


%% WAVE

load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/grid bias testing/circletest/wave_Kc1_5e-08.mat')

% ANALYTICAL
Area0 = 25*pi;
r0 = 5;
k_wave = p.Kcoast;
t = 1;
Areat_wave = pi*((r0/(1-r0*k_wave*p.dx*3/4*pi*t)).^2);
dArea_wave_analytical = Areat_wave - Area0;

% MODELED
wave_init = ~(output(:,:,1));
wave_final = g.U<1;
wave_init_area = length(find(wave_init))./((length(wave_init)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
wave_final_area = length(find(wave_final))./((length(wave_final)).^2)*(p.dx*p.Nx*p.dy*p.Ny);
dArea_wave_model = wave_final_area - wave_init_area;