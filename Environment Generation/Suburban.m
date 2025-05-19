% Parameters for Suburban Environment with Gradual Altitude Variations
terrainSize = 200; % Size of the terrain (e.g., 100x100)
numBuildings = 120; % Fewer buildings for suburban areas
minBuildingSize = 10; % Minimum size of a suburban house
maxBuildingSize = 10; % Maximum size of a suburban house
maxBuildingHeight = 20; % Maximum height for suburban buildings (shorter than city buildings)
altitudeScale = 10; % Increased scale for more spread-out altitude variations

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

% Randomly place suburban houses on the terrain
for i = 1:numBuildings
    % Random building size
    buildingWidth = randi([minBuildingSize, maxBuildingSize]);
    buildingDepth = randi([minBuildingSize, maxBuildingSize]);
    buildingHeight = randi([5, maxBuildingHeight]); % Lower height range for houses
    
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
terrainStruct.end = [200; 200; 40];
terrainStruct.n = 7;
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

save('terrainStruct.mat', 'terrainStruct');

% Plot the terrain with suburban houses
figure;
surf(X, Y, baseTerrain);
xlabel('X');
ylabel('Y');
zlabel('Z');
%title('Suburban Terrain with Gradual Altitude Variations');
%colorbar;
%shading interp; % Optional: for smoother shading

% Display the struct
disp(terrainStruct);
