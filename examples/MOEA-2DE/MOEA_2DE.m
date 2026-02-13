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
    % Problem Definition
    Gen=500;          % Maximum Number of Iterations
    pop=20;           % Population Size
    pm=0.5;           % probability of mutation
    pc=0.9;           % probability of crossover
    M=2;              % objective number
    CR=0.3;F=0.3;     % DE parameter
      
    Runs = 30;         
       
    j=5;
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
    for fld = 1:Runs    
        if not (isfolder(sprintf('Population%s/Run_%s',tempp,num2str(fld))))
            mkdir(sprintf('Population%s/Run_%s',tempp,num2str(fld)));
        end
    end
    bestScores = [];
    for run = 1 : Runs
        disp(run)
        eval = Gen * pop;
        pro = 0.3;
        adjustGeneration = 1;
        gen_hv = [];
        % population Initialization
        for i = 1 : pop
            population(i) = Chromosome(model);
            child(i) = Chromosome(model);
            population(i) = initialize(population(i),model);
            population(i) = evaluate(population(i));
        end
        [~,FrontNo,CrowdDis] = EnvironmentalSelection(population,pop,pop,M);
        dim = size(population(1).rnvec,1); % Number of waypoints
        Score(1,1) = 0;
        Score(1,2) = 0;
        
        strongDim = zeros(1,dim);
        b = strongDim;
        parent = Chromosome(model);
        parent.rnvec(1,:) = model.start;
        parent.rnvec(model.n,:) = model.end;
        parent.rnvec(2:dim-1,1) = [2:dim-1]*model.xmax/model.n;
        parent.rnvec(2:dim-1,2) = [2:dim-1]*model.xmax/model.n;
        parent.path = testBspline([parent.rnvec(:,1)';parent.rnvec(:,2)'],model.xmax)';
        parent = adjust_constraint_turning_angle(parent,model);
        parent = evaluate(parent);
        %key dimension exploration
        times = 3;
        varRange = times*model.xmax/model.n;
        for j = 1 : dim - 2
            testIndiv = Chromosome(model);
            testIndiv.rnvec = parent.rnvec;
            rnd = randi(2);
            testIndiv.rnvec(j+1,rnd) = parent.rnvec(j+1,rnd) + varRange; 
            testIndiv.rnvec = Evolve.check_boundary(testIndiv.rnvec,j+1,model);
            % check constraint
            testIndiv.path = testBspline([testIndiv.rnvec(:,1)';testIndiv.rnvec(:,2)'],model.xmax)';
            testIndiv = adjust_constraint_turning_angle(testIndiv,model);
            testIndiv = evaluate(testIndiv);
            dimDiff(j,1) = abs(parent.objs(1) - testIndiv.objs(1));
            dimDiff(j,2) = abs(parent.objs(2) - testIndiv.objs(2));
        end
        v = mean(dimDiff); % Objective values
        a = find(dimDiff(:,1)>v(1));
        b(a)=1;
        a = find(dimDiff(:,2)>v(2));
        b(a)=1;
        strongDim = b;
        strongDim(1)=0;
        strongDim(end)=0;
        eval = eval - (dim-1);
        
        % Main Loop
        i = 0;
        while 1
            lastpopulation = population;
            eliteCluster = population(FrontNo == 1);
            % Create offspring population
            MatingPool = TournamentSelection(2,pop,FrontNo,-CrowdDis);
            for j = 1 : 2 : pop-1 % pop: Population size
                if rand(1)<pro 
                    parent1 = eliteCluster(randi(length(eliteCluster)));
                    [child(j)] = Evolve.dimExplore(parent1,dim,model,strongDim,F);
                    parent2 = eliteCluster(randi(length(eliteCluster)));
                    [child(j+1)] = Evolve.dimExplore(parent2,dim,model,strongDim,F);
                else 
                    [child(j),child(j+1)] = Evolve.binary_crossover(population(MatingPool(j)),population(MatingPool(j+1)),nVar-1,pc,randi(2));
                    if rand < pm                        
                        child(j) = Evolve.mutation(child(j),nVar,model);
                    end
                    if rand < pm                        
                        child(j+1) = Evolve.mutation(child(j+1),nVar,model);                        
                    end
                end
                child(j).rnvec = sortrows(child(j).rnvec,randi(2));
                child(j+1).rnvec = sortrows(child(j+1).rnvec,randi(2));
                child(j).path = testBspline([child(j).rnvec(:,1)';child(j).rnvec(:,2)'],model.xmax)';
                child(j) = adjust_constraint_turning_angle(child(j),model);
                child(j) = evaluate(child(j));
                child(j+1).path = testBspline([child(j+1).rnvec(:,1)';child(j+1).rnvec(:,2)'],model.xmax)';
                child(j+1) = adjust_constraint_turning_angle(child(j+1),model);
                child(j+1) = evaluate(child(j+1));
        
                i = i + 2;
                if i >= eval
                    break;
                end
            end
           
            [population,FrontNo,CrowdDis] = EnvironmentalSelection([population,child],pop,pop*2,M);
            
            obj = [population.objs];
            PopObj = reshape(obj,M,length(population))';
            g = ceil(i/pop + 1); 
            difference(g) = KLD(dim,lastpopulation,population);
            if g > 1
                if mod(g-1,adjustGeneration) == 0
                    j = g-adjustGeneration;
                    upGap = 0;
                    downGap = 0;
                    while j < g
                        if difference(j+1) - difference(j) > 0
                            upGap = upGap + difference(j+1) - difference(j);
                        else
                            downGap = downGap + difference(j) - difference(j+1);
                        end
                        j = j + 1;
                    end
                    if upGap > downGap
                        pro = 1;
                    else
                        pro = 0;
                    end
                end
            end
            
            obj = [population.objs];
            PopObj = reshape(obj,M,length(population))';
            [cur_hv] = [calMetirc(1,PopObj,problemIndex),calMetirc(2,PopObj,problemIndex)];
            gen_hv = [gen_hv,cur_hv];
            
            if i >= eval
                break;
            end
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
