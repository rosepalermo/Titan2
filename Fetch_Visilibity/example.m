% Visibility algorithm/find fetch demo
% Rose Palermo - Updated May 5 2021
% This demo shows a calculation of the fetch around a lake using the
% visilibity1 algorithm, assuming isotropic winds.

% INPUT: shoreline - shoreline of a lake, formatted such that it is one cell with the
% first element the shoreline points in order and subsequent elements are
% ordered island points. Shoreline ordering should be CW and island
% ordering should be CCW.

% OUTPUT: All fetch and associated cos(phi-theta) for each point on the
% shoreline. 

% (Una, you'll want to use the cos(phi-theta) to find phi and theta for
% limiting wind directions
%% before running the demo, need to mex things in the src folder (see ReadMe)
% mex -setup
% mex -v in_environment.cpp visilibity.o
% mex -v visibility_polygon.cpp visilibity.o
%% load example input
load('example_input.mat')

%% run

% parameters for visibility functions
eps = 0.1;
epsilon = 1e-4;
snap_distance = 0.05; % Looks like we get some invalid (empty array) visibility polygons if snap_distance >= eps.

[WaveArea, Fetch_dist, cosang, V] = fetch_VisiLibity(shorelines,eps,epsilon,snap_distance);


% plot total wave area
figure()
for k = 1:length(shorelines)
scatter3(shorelines{k}(:,1),shorelines{k}(:,2),WaveArea{k}(:,1),[],WaveArea{k}(:,1),'.')
hold on
view(2)
end

figure()
plot_env(shorelines, [shorelines{1}(1,1) shorelines{1}(1,2)],V);
