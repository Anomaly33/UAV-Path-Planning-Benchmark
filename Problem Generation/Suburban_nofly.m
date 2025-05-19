% Parameters for Suburban Environment with Gradual Altitude Variations
terrainSize = 200; % Size of the terrain (e.g., 100x100)
numBuildings = 110; % Fewer buildings for suburban areas
minBuildingSize = 10; % Minimum size of a suburban house
maxBuildingSize = 10; % Maximum size of a suburban house
maxBuildingHeight = 20; % Maximum height for suburban buildings (shorter than city buildings)
altitudeScale = 10; % Increased scale for more spread-out altitude variations
noFlyZoneRadius = 20; % Radius of the no-fly zone (cylinder)
noFlyZoneHeight = 60; % Height of the no-fly zone (higher than buildings)

% Generate smooth gradual altitude variations (2-3 significant hills)
[X, Y] = meshgrid(1:terrainSize, 1:terrainSize);

% Use a combination of sinusoidal functions with larger wavelengths for smooth, gradual hills
baseTerrain = altitudeScale * ( ...
    sin(0.005 * X) + sin(0.005 * Y) + ... % Large smooth hills
    0.5 * sin(0.03 * X) .* sin(0.03 * Y) ... % Additional undulations
);

% Normalize the terrain to a smooth gradual range
baseTerrain = (baseTerrain - min(baseTerrain(:))); % Shift the terrain so the minimum value is 0
baseTerrain = baseTerrain / max(baseTerrain(:)); % Normalize to a range of 0 to 1
baseTerrain = baseTerrain * maxBuildingHeight * 2; % Scale for significant gradual altitude variation

% Define the no-fly zone as a cylinder in the terrain
noFlyZoneCenterX = terrainSize / 2;
noFlyZoneCenterY = terrainSize / 2;
[X_grid, Y_grid] = meshgrid(1:terrainSize, 1:terrainSize);
distanceFromCenter = sqrt((X_grid - noFlyZoneCenterX).^2 + (Y_grid - noFlyZoneCenterY).^2);

% Add the cylindrical no-fly zone to the base terrain
noFlyZoneMask = distanceFromCenter <= noFlyZoneRadius; % Create a mask for the cylinder
baseTerrain(noFlyZoneMask) = noFlyZoneHeight; % Set the height of the no-fly zone

% Randomly place suburban houses on the terrain
for i = 1:numBuildings
    % Random building size
    buildingWidth = randi([minBuildingSize, maxBuildingSize]);
    buildingDepth = randi([minBuildingSize, maxBuildingSize]);
    buildingHeight = randi([5, maxBuildingHeight]); % Lower height range for houses
    
    % Random position (ensure the building fits within the terrain)
    % startX = randi([1, terrainSize - buildingWidth + 1]);
    % startY = randi([1, terrainSize - buildingDepth + 1]);
    
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

% Create the struct
terrainStruct = struct();
terrainStruct.start = [1; 1; 1];
terrainStruct.end = [200; 200; 40];
terrainStruct.n = 10;
terrainStruct.xmin = 1;
terrainStruct.xmax = terrainSize;
terrainStruct.ymin = 1;
terrainStruct.ymax = 200;
terrainStruct.zmin = 0; % Ensure the minimum z-value is set to ground level
terrainStruct.zmax = 50; % Adjust zmax based on the generated terrain
terrainStruct.X = X(1,:);
terrainStruct.Y = Y(:,1);
terrainStruct.H = baseTerrain;
terrainStruct.safeH = 5; % Adjusted safe height for suburban environment
terrainStruct.theta = 30; % Adjust as needed
terrainStruct.nofly_c = [noFlyZoneCenterX,noFlyZoneCenterY];
terrainStruct.nofly_r = noFlyZoneRadius;
terrainStruct.nofly_h = noFlyZoneHeight;

save('terrainStruct.mat', 'terrainStruct');

% Plot the terrain with suburban houses
figure;
surf(X, Y, baseTerrain);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Suburban Terrain with Gradual Altitude Variations');
colorbar;
%shading interp; % Optional: for smoother shading
hold on

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
% Display the struct
disp(terrainStruct);
