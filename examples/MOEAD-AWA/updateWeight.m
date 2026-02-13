function [Population,W] = updateWeight(Population,W,Z,EP,nus)
% Delete overcrowded subproblems and add new subproblems

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    obj = [Population.objs];
    PopObj = reshape(obj,2,length(Population))';
    [N,M] = size(PopObj);

    %% Update the current population by EP
    % Calculate the function value of each solution in Population or EP on
    % each subproblem in W
    Combine = [Population,EP];
    com_obj = [Combine.objs];
    ComObj = reshape(com_obj,2,length(Combine))';
    CombineObj = abs(ComObj-repmat(Z,length(Combine),1));
    g = zeros(length(Combine),size(W,1));
    for i = 1 : size(W,1)
        g(:,i) = max(CombineObj.*repmat(W(i,:),length(Combine),1),[],2);
    end
    % Choose the best solution for each subproblem
    [~,best]   = min(g,[],1);
    Population = Combine(best);
    
    %% Delete the overcrowded subproblems
    Dis = pdist2(PopObj,PopObj);
    Dis(logical(eye(length(Dis)))) = inf;
    Del = false(1,length(Population));
    while sum(Del) < min(nus,length(EP))
        Remain = find(~Del);
        subDis = sort(Dis(Remain,Remain),2);
        [~,worst] = min(prod(subDis(:,1:min(M,length(Remain))),2));
        Del(Remain(worst)) = true;
    end
    Population = Population(~Del);
    W = W(~Del,:);
    
    %% Add new subproblems
    % Determine the new solutions be added
    Combine  = [Population,EP];
    Selected = false(1,length(Combine));
    Selected(1:length(Population)) = true;
    com_obj = [Combine.objs];
    ComObj = reshape(com_obj,2,length(Combine))';
    Dis = pdist2(ComObj,ComObj);
    Dis(logical(eye(length(Dis)))) = inf;
    while sum(Selected) < min(N,length(Selected))
        subDis = sort(Dis(~Selected,Selected),2);
        [~,best] = max(prod(subDis(:,1:min(M,size(subDis,2))),2));
        Remain = find(~Selected);
        Selected(Remain(best)) = true;
    end
    % Add new subproblems
    ep_obj = [EP(Selected(length(Population)+1:end)).objs];
    ep_obj = reshape(ep_obj,2,length(EP(Selected(length(Population)+1:end))))';
    % newObjs = EP(Selected(length(Population)+1:end)).objs;
    newObjs = ep_obj;
    temp    = 1./(newObjs-repmat(Z,size(newObjs,1),1));
    W = [W;temp./repmat(sum(temp,2),1,size(temp,2))];
    % Add new solutions
    Population = Combine(Selected);
end