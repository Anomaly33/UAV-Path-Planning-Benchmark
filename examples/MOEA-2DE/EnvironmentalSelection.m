function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N,pop,M)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    obj = [Population.objs];
    obj = reshape(obj,M,pop)';
    conss = [Population.cons]';
    
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(obj,conss,N);

    % x = obj(:,1);
    % y = obj(:,2);
    % % Create the scatter plot
    % figure;
    % scatter(x, y,70, 'filled');
    % 
    % % Annotate each point with its index
    % numPoints = length(x);
    % hold on;
    % for i = 1:numPoints
    %     text(x(i), y(i), num2str(FrontNo(i)), 'FontSize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    % end
    % hold off;
    % 
    % % Label the axes
    % xlabel('X-axis');
    % ylabel('Y-axis');
    % title('Scatter Plot with Numbered Points');
    % 
    % % Generate a unique filename
    % filename = sprintf('Plots/im.png');
    % %filepath = fullfile(savePath, filename);
    % 
    % % Save the figure
    % saveas(gcf, filename);
    % 
    % % Close the figure to free up memory
    % close(gcf);

    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(obj,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end