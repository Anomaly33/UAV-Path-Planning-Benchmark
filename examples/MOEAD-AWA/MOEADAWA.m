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
    boundary = [MaxValue;MinValue];
    
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
        Score(1,1) = 0;
        Score(1,2) = 0;
        FunctionValue=[];
        gen_hv = [];
        %% Generate the weight vectors
        rate_update_weight = 0.05;
        rate_evol = 0.8;
        wag = 100;
        [W,pop] = UniformPoint(pop,M);
        W = 1./W./repmat(sum(1./W,2),1,size(W,2));
        nr = ceil(pop/100);
        nEP = ceil(pop*1.5);
        T = ceil(pop/10);

        %% Detect the neighbours of each solution
        B = pdist2(W,W);
        [~,B] = sort(B,2);
        B = B(:,1:T);
        %% Generate random population
        for i = 1 : pop
            population(i) = Chromosome(model);
            population(i) = initialize(population(i),model);
            population(i) = evaluate(population(i));
        end
        objs = [population.objs];
        obj = reshape(objs,M,length(population))';
        Z = min(obj,[],1);
        Pi = ones(pop,1);
        oldObj = max(abs((obj-repmat(Z,pop,1)).*W),[],2);
        %% Optimization
        EP = [];
        gen = 0;
        
        while gen<Generations
            if ~mod(ceil(gen*pop/pop),10)
                objs = [population.objs];
                obj = reshape(objs,M,length(population))';
                newObj    = max(abs((obj-repmat(Z,pop,1)).*W),[],2);
                DELTA     = (oldObj-newObj)./oldObj;
                Temp      = DELTA <= 0.001;
                Pi(~Temp) = 1;
                Pi(Temp)  = (0.95+0.05*DELTA(Temp)/0.001).*Pi(Temp);
                oldObj    = newObj;
            end
            for subgeneration = 1 : 5
                Bounday = find(sum(W<1e-3,2)==1)';
                I = [Bounday,TournamentSelection(10,floor(pop/5)-length(Bounday),-Pi)];
                for i = 1 : length(I)
                    if rand < 1%0.9
                        P = B(I(i),randperm(size(B,2)));
                    else
                        P = randperm(pop);
                    end
                    offspring(i) = F_operator(population(i),population(P(1:2)),boundary,model);
                    Z = min(Z,offspring(i).objs);

                    pop_obj = [population(P).objs];
                    pop_obj = reshape(pop_obj,M,length(P))';
                    g_old = max(abs(pop_obj-repmat(Z,T,1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(offspring(i).objs-Z),T,1).*W(P,:),[],2);
                    population(P(find(g_old>=g_new,nr))) = offspring(i);
                end
            end
            if gen*pop >=rate_evol*Generations*pop
            % Adaptive weight adjustment
                if isempty(EP)
                    EP = updateEP(population,offspring,nEP);
                else
                    EP = updateEP(EP,offspring,nEP);
                end
                if ~mod(ceil(gen*pop/pop),wag/5)
                    [population,W] = updateWeight(population,W,Z,EP,rate_update_weight*pop);
                end
            end
            obj = [population.objs];
            PopObj = reshape(obj,M,length(population))';
            [cur_hv] = [calMetirc(1,PopObj,problemIndex),calMetirc(2,PopObj,problemIndex)];
            gen_hv = [gen_hv;cur_hv];
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
                dt_sv.path = population(1,i).path;
                dt_sv.objs = population(1,i).objs;
                save(sprintf('Population%s/Run_%s/bp_%s.mat',tempp,num2str(run),num2str(pp)),'dt_sv');
                pp=pp+1;
        end
    end
    save(sprintf('Population%s/final_hv.mat',tempp),'bestScores');
    clearvars -except files pro_num
end