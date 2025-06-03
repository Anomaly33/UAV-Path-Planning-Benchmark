clear all; clc;

d = dir(fullfile('problems','*.mat'));
temp = [d.name];
files = split(temp,'.mat');
files = files(1:end-1);

for pro_num = 1:numel(files)
    pro_nme = strcat(files(pro_num),'.mat');
    pro_nme = pro_nme{:};
    data = load(pro_nme);
    tempp = split(pro_nme,'terrainStruct');
    tempp = tempp{2};
    tempp = split(tempp,'.mat');
    tempp = tempp{1};
    fprintf('Current Problem: %s',tempp(2:end));
    
    terrainSize = 200;

    % data = load('terrainStruct_s_110_20_nofly');
    terrain = data.terrainStruct.H;
    
    X = data.terrainStruct.X;
    Y = data.terrainStruct.Y;
    
    % Visualize the terrain using the surf function
    fig = figure('Visible', 'off'); % Create a figure handle
    surf(X, Y, terrain);
    
    % Set view orientation
    view(-10, 30);
    
    % Add labels and a colorbar
    xlabel('X');
    ylabel('Y');
    zlabel('Height (Z)');
    cb = colorbar; % Get the colorbar handle
    cb.Title.String = 'Altitude'; % Set the colorbar title
    cb.Title.Position = [0, 388.5, 0];
    cb.Title.VerticalAlignment = 'bottom';
    % set(gca, 'Position', [0, 0, 0, 0.1]);
    % Optional: Uncomment for smoother shading
    % shading interp;
    
    % Tighten the layout
    axis tight;
    set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02)); % Minimize margins
    
    % Set figure properties for high-quality saving
    set(fig, 'PaperPositionMode', 'auto'); % Match layout to screen
    set(fig, 'Units', 'inches', 'Position', [0, 0, 10, 6]); % Adjust size (8x6 inches)
    
    % Save the figure as a high-quality PNG
    print(fig, sprintf('figs/%s.png',tempp), '-dpng', '-r300'); % '-r300' for 300 DPI
    close(fig)
end
