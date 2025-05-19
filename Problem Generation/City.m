% Parameters
terrainSize = 200; % Size of the terrain (e.g., 100x100)
numBuildings = 120; % Number of buildings
minBuildingSize = 3; % Minimum size of a building
maxBuildingSize = 15; % Maximum size of a building
maxBuildingHeight = 305;%3500; % Maximum height of a building
altitudeScale = 1; % Scale for the terrain altitude variations for subtle changes

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

% Randomly place buildings on top of the base terrain
for i = 1:numBuildings
    % Random building size
    buildingWidth = randi([minBuildingSize, maxBuildingSize]);
    buildingDepth = randi([minBuildingSize, maxBuildingSize]);
    buildingHeight = randi([1, maxBuildingHeight]);
    
    % Random position (ensure the building fits within the terrain)
    startX = randi([1, terrainSize - buildingWidth + 1]);
    startY = randi([1, terrainSize - buildingDepth + 1]);
    
    % Calculate the average height of the terrain where the building will be placed
    avgHeight = mean(mean(baseTerrain(startY:startY+buildingDepth-1, startX:startX+buildingWidth-1)));
    
    % Place the building on the terrain by adding height
    baseTerrain(startY:startY+buildingDepth-1, startX:startX+buildingWidth-1) = avgHeight + buildingHeight;
end

% Create the struct
terrainStruct = struct();
terrainStruct.start = [1; 1; 1];
terrainStruct.end = [200; 200; 100];
terrainStruct.n = 7;
terrainStruct.xmin = 1;
terrainStruct.xmax = terrainSize;
terrainStruct.ymin = 1;
terrainStruct.ymax = 200;
terrainStruct.zmin = 1;%min(baseTerrain(:));
terrainStruct.zmax = 120;%3500; % Adjust this value based on your data
terrainStruct.X = X(1,:);
terrainStruct.Y = Y(:,1);
terrainStruct.H = baseTerrain;
terrainStruct.safeH = 90;%maxBuildingHeight;
terrainStruct.theta = 30; % Adjust this value based on your needs

save('terrainStruct.mat', 'terrainStruct');

% Plot the terrain with buildings
figure;
surf(X, Y, baseTerrain);
xlabel('X');
ylabel('Y');
zlabel('Z');
%title('Subtle Terrain with Buildings and Natural Hills');
%colorbar;
%shading interp; % Optional: for smoother shading

% Display the struct
disp(terrainStruct);

