clear all;clc;format compact;tic;
%% Automatically load the problem files and run the algo iteratively
d = dir(fullfile('problems','*.mat'));
temp = [d.name];
files = split(temp,'.mat');
files = files(1:end-1);

for pro_num = 1:numel(files)
    pro_nme = strcat(files(pro_num),'.mat');
    pro_nme = pro_nme{:};
    load(pro_nme);
    tempp = split(pro_nme,'terrainStruct');
    tempp = tempp{2};
    tempp = split(tempp,'.mat');
    tempp = tempp{1};
    fprintf('Current Problem: %s',tempp(2:end));
    j=5;    %Choosing the terrain
    if j == 1
        model = peak1;        
        problemIndex = 1;
    elseif j == 2
        model = peak2;
        problemIndex = 1;
    elseif j == 3
        model = Rugged1;
        problemIndex = 2;
    elseif j ==4
        model = Rugged2;
        problemIndex = 2;
    else
        model = terrainStruct;
        problemIndex = 3;
    end
    nVar=model.n;       
    
    MinValue   = [model.xmin,model.ymin,model.zmin];%zeros(1,D);
    MaxValue   = [model.xmax,model.ymax,model.zmax];%ones(1,D);
    Boundary = [MaxValue;MinValue];
    
    Generations = 500;
    pop = 20; 
    % N=3;
    M=2;
    Runs=30;
    bestScores=[];
    type = 4;
    
    for fld = 1:Runs    
        if not (isfolder(sprintf('Population%s/Run_%s',tempp,num2str(fld))))
            mkdir(sprintf('Population%s/Run_%s',tempp,num2str(fld)));
        end
    end
    for run=1:Runs
        disp(run)
        g=2;
        Score(1,1)=0;
        Score(1,2)=0;
        FunctionValue=[];
        gen_hv = [];
        %% Generate the weight vectors
        [W,pop] = UniformPoint(pop,M);
        T = ceil(pop/10);

        %% Detect the neighbours of each solution
        B = pdist2(W,W);
        [~,B] = sort(B,2);
        B = B(:,1:T);
        %% Generate random population
        for i = 1 : pop
            population(i) = Chromosome(model);% Set up Parent population [1 1 300;0 0 0;....]
            % child(i) = Chromosome(model);% Set up Child Population
            population(i) = initialize(population(i),model); % Initializes all the solutions
            population(i) = evaluate(population(i));
        end
        objs = [population.objs];
        obj = reshape(objs,M,length(population))';
        Z = min(obj,[],1);
        %% Optimization
        gen = 0;
        while gen<Generations
            % disp(gen)
            for i = 1:pop
                P = B(i,randperm(size(B,2)));
                offspring = F_operator(population(i),population(P(1:2)),Boundary,model);
                Z = min(Z,offspring.objs);
                pop_obj = [population(P).objs];
                pop_obj = reshape(pop_obj,M,length(P))';
                if type==1
                    normW   = sqrt(sum(W(P,:).^2,2));
                    normP   = sqrt(sum((pop_obj-repmat(Z,T,1)).^2,2));
                    normO   = sqrt(sum((offspring.objs-Z).^2,2));
                    CosineP = sum((pop_obj-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                    CosineO = sum(repmat(offspring.objs-Z,T,1).*W(P,:),2)./normW./normO;
                    g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                    g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                    population(P(g_old>=g_new)) = offspring;
                elseif type==2
                    g_old = max(abs(pop_obj-repmat(Z,T,1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(offspring.objs-Z),T,1).*W(P,:),[],2);
                elseif type==3
                    Zmax  = max(pop_obj,[],1);
                    g_old = max(abs(pop_obj-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                    g_new = max(repmat(abs(offspring.objs-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
                else
                    g_old = max(abs(pop_obj-repmat(Z,T,1))./W(P,:),[],2);
                    g_new = max(repmat(abs(offspring.objs-Z),T,1)./W(P,:),[],2);
                end
                population(P(g_old>=g_new)) = offspring;
            end
            obj = [population.objs];
            PopObj = reshape(obj,M,length(population))';
            [cur_hv] = [calMetirc(1,PopObj,problemIndex),calMetirc(2,PopObj,problemIndex)];
            gen_hv = [gen_hv,cur_hv];
            gen=gen+1;
        end
        save(sprintf('Population%s/Run_%s/gen_hv.mat',tempp,num2str(run)),'gen_hv');
        obj = [population.objs];
        PopObj = reshape(obj,M,length(population))';

        [Score(1,1)] = calMetirc(1,PopObj,problemIndex);
        [Score(1,2)] = calMetirc(2,PopObj,problemIndex);
        bestScores = [bestScores,Score];
        pp=1;
        for i = 1:size(population,2)
            %filen = sprintf('Population/Run_%s/bp_%s.mat',num2str(run),num2str(i));
            % if ~PopObj(i,3)>0
                dt_sv.path = population(1,i).path;
                dt_sv.objs = population(1,i).objs;
                save(sprintf('Population%s/Run_%s/bp_%s.mat',tempp,num2str(run),num2str(pp)),'dt_sv');
                pp=pp+1;
            % end
        end
    end
    save(sprintf('Population%s/final_hv.mat',tempp),'bestScores');
    clearvars -except files pro_num
end