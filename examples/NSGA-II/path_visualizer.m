clear all;clc;
% Define the terrain size
terrainSize = 200;

data = load('terrainStruct_m_200_25nc_nofly.mat');

run = 1;

d = dir([sprintf('Population_m_200_25nc_nofly/Run_%s/',num2str(run)), '/*.mat' ]);
length(d)
sid = randi([1,length(d)-1],1); % Randomly selected individual from pareto front of a selected run

temp = load(sprintf('Population_m_200_25nc_nofly/Run_%s/bp_%s.mat',num2str(run),num2str(sid)));
sol = temp.dt_sv.path;

% Create a 200x200 matrix representing the terrain heights (example data)
% You should replace this with your actual 200x200 matrix
terrain = data.terrainStruct.H; % Random heights for demonstration

% Create meshgrid for X and Y coordinates
% [X, Y] = meshgrid(1:terrainSize, 1:terrainSize);

X = data.terrainStruct.X;
Y = data.terrainStruct.Y;

% Visualize the terrain using the surf function
figure;
surf(X, Y, terrain);
xlabel('X');
ylabel('Y');
zlabel('Height (Z)');
title('S-180');
colorbar;
%shading interp; % Optional: for smoother shading

hold on;
uav_x = sol(:,1)';
uav_y = sol(:,2)';
uav_z = sol(:,3)';

plot3(uav_x, uav_y, uav_z, 'r-', 'LineWidth', 2);

legend('Terrain', 'UAV Path');

% Release the hold
hold off;
