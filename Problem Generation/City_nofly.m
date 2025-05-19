% Parameters
terrainSize = 200; % Size of the terrain
numBuildings = 100; % Number of buildings
minBuildingSize = 3; % Minimum size of a building
maxBuildingSize = 15; % Maximum size of a building
maxBuildingHeight = 305; % Maximum height of a building
altitudeScale = 1; % Scale for the terrain altitude variations
noFlyZoneRadius = 20; % Radius of the no-fly zone (cylinder)
noFlyZoneHeight = 400; % Height of the no-fly zone (higher than buildings)

% Generate the base terrain with subtle hills
[X, Y] = meshgrid(1:terrainSize, 1:terrainSize);
baseTerrain = altitudeScale * ( ...
    sin(0.1 * X) .* sin(0.1 * Y) + ...
    0.5 * sin(0.2 * X) .* sin(0.2 * Y) + ...
    0.25 * sin(0.3 * X) .* sin(0.3 * Y) + ...
    0.125 * sin(0.4 * X) .* sin(0.4 * Y) ...
);

% Normalize the terrain to a specific range
baseTerrain = (baseTerrain - min(baseTerrain(:))) / (max(baseTerrain(:)) - min(baseTerrain(:)));
baseTerrain = baseTerrain * maxBuildingHeight / 10; % Scale to a smaller range for subtle variations

% Define the no-fly zone as a cylinder in the terrain
noFlyZoneCenterX = terrainSize / 2;
noFlyZoneCenterY = terrainSize / 2;
[X_grid, Y_grid] = meshgrid(1:terrainSize, 1:terrainSize);
distanceFromCenter = sqrt((X_grid - noFlyZoneCenterX).^2 + (Y_grid - noFlyZoneCenterY).^2);

% Add the cylindrical no-fly zone to the base terrain
noFlyZoneMask = distanceFromCenter <= noFlyZoneRadius; % Create a mask for the cylinder
baseTerrain(noFlyZoneMask) = noFlyZoneHeight; % Set the height of the no-fly zone

% Randomly place buildings on top of the base terrain (excluding the no-fly zone)
for i = 1:numBuildings
    % Random building size
    buildingWidth = randi([minBuildingSize, maxBuildingSize]);
    buildingDepth = randi([minBuildingSize, maxBuildingSize]);
    buildingHeight = randi([1, maxBuildingHeight]);
    
    % Random position (ensure the building fits within the terrain and not in no-fly zone)
    validPosition = false;
    while ~validPosition
        startX = randi([1, terrainSize - buildingWidth + 1]);
        startY = randi([1, terrainSize - buildingDepth + 1]);
        if ~any(noFlyZoneMask(startY:startY+buildingDepth-1, startX:startX+buildingWidth-1), 'all')
            validPosition = true;
        end
    end
    
    % Calculate the average height of the terrain where the building will be placed
    avgHeight = mean(mean(baseTerrain(startY:startY+buildingDepth-1, startX:startX+buildingWidth-1)));
    
    % Place the building on the terrain by adding height
    baseTerrain(startY:startY+buildingDepth-1, startX:startX+buildingWidth-1) = avgHeight + buildingHeight;
end

% Create the struct with no-fly zone included
terrainStruct = struct();
terrainStruct.start = [1; 1; 1];
terrainStruct.end = [200; 200; 100];
terrainStruct.n = 7;
terrainStruct.xmin = 1;
terrainStruct.xmax = terrainSize;
terrainStruct.ymin = 1;
terrainStruct.ymax = 200;
terrainStruct.zmin = 1;
terrainStruct.zmax = 120;
terrainStruct.X = X(1,:);
terrainStruct.Y = Y(:,1);
terrainStruct.H = baseTerrain;
terrainStruct.safeH = 90;
terrainStruct.theta = 30;
terrainStruct.nofly_c = [noFlyZoneCenterX,noFlyZoneCenterY];
terrainStruct.nofly_r = noFlyZoneRadius;
terrainStruct.nofly_h = noFlyZoneHeight;

save('terrainStruct.mat', 'terrainStruct');

% Plot the terrain with buildings and no-fly zone as a mesh
figure;
surf(X, Y, baseTerrain);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Terrain with Buildings and No-Fly Zone');
colorbar;
shading interp; % Optional: for smoother shading
hold on;

% Plot the no-fly zone as a transparent mesh (optional for visualization)
theta = linspace(0, 2*pi, 100);
circleX = noFlyZoneCenterX + noFlyZoneRadius * cos(theta);
circleY = noFlyZoneCenterY + noFlyZoneRadius * sin(theta);
zMesh = linspace(0, noFlyZoneHeight, 50);
[CircleX, ZMesh] = meshgrid(circleX, zMesh);
[CircleY, ~] = meshgrid(circleY, zMesh);

mesh(CircleX, CircleY, repmat(zMesh', 1, length(theta)), 'FaceAlpha', 0.4, 'EdgeColor', 'r'); % Semi-transparent

% Adjust plot limits
xlim([1, terrainSize]);
ylim([1, terrainSize]);
zlim([0, noFlyZoneHeight]);

hold off;
