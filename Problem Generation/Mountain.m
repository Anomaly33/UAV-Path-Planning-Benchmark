% Parameters for Isolated Mountainous Terrain
terrainSize = 200; % Size of the terrain (e.g., 200x200)
numPeaks = 200; % Number of isolated mountain peaks
peakHeight = 300; % Maximum height of the mountain peaks
valleyDepth = 0; % Depth of the valleys
altitudeScale = 8; % Scale for the terrain altitude variations

% Generate the base terrain with random noise for natural ruggedness
[X, Y] = meshgrid(1:terrainSize, 1:terrainSize);
baseTerrain = altitudeScale * ( ...
    sin(0.03 * X) .* sin(0.03 * Y) + ... % Low frequency for large-scale undulations
    0.5 * sin(0.07 * X) .* sin(0.07 * Y) + ... % Medium frequency variations
    0.25 * sin(0.1 * X) .* sin(0.1 * Y) + ... % Small-scale undulations
    0.1 * rand(terrainSize) ... % Random noise for rugged terrain
);

% Add isolated mountain peaks at random locations
for i = 1:numPeaks
    peakX = randi([17, terrainSize - 17]); % Random X position (avoid edges)
    peakY = randi([17, terrainSize - 17]); % Random Y position (avoid edges)
    peakRadius = randi([15, 30]); % Random radius for each peak

    % Create circular peaks with a Gaussian shape
    for x = max(1, peakX - peakRadius):min(terrainSize, peakX + peakRadius)
        for y = max(1, peakY - peakRadius):min(terrainSize, peakY + peakRadius)
            distanceToPeak = sqrt((x - peakX)^2 + (y - peakY)^2);
            if distanceToPeak <= peakRadius
                baseTerrain(y, x) = baseTerrain(y, x) + ...
                    peakHeight * exp(-distanceToPeak^2 / (2 * (peakRadius / 2)^2));
            end
        end
    end
end

% Normalize the terrain to a specific range between valleys and peaks
baseTerrain = (baseTerrain - min(baseTerrain(:))) / (max(baseTerrain(:)) - min(baseTerrain(:)));
baseTerrain = baseTerrain * (peakHeight - valleyDepth) + valleyDepth;

% Create the struct for terrain data
terrainStruct = struct();
terrainStruct.start = [1; 1; 1]; % Starting position in the terrain
terrainStruct.end = [200; 200; 1]; % Goal at a higher altitude
terrainStruct.n = 7;
terrainStruct.xmin = 1;
terrainStruct.xmax = terrainSize;
terrainStruct.ymin = 1;
terrainStruct.ymax = terrainSize;
terrainStruct.zmin = valleyDepth;
terrainStruct.zmax = peakHeight;
terrainStruct.X = X(1,:);
terrainStruct.Y = Y(:,1);
terrainStruct.H = baseTerrain;
terrainStruct.safeH = 10; % Safe height for flying
terrainStruct.theta = 30; % Adjusted for steeper slopes in a mountainous terrain

save('terrainStruct.mat', 'terrainStruct');

% Plot the terrain with spread-apart mountains
figure;
surf(X, Y, baseTerrain);
xlabel('X');
ylabel('Y');
zlabel('Z');
%title('Isolated Mountain Terrain with Spread-Apart Peaks');
%colorbar;
%shading interp; % Optional: for smoother shading

% Display the struct
disp(terrainStruct);
