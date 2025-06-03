function newpop = F_operator(population,MatingPool,Boundary,model)
% This function generates a new population by genetic operators
    
    N = size(MatingPool,1);
    D = 3;
%-----------------------------------------------------------------------------------------
% Parameters setting
    ProC = 1;       % The probability of crossover
    ProM = 1/D;     % The probability of mutation
    DisC = 20;   	% The parameter of crossover
    DisM = 20;   	% The parameter of mutation
%-----------------------------------------------------------------------------------------
% Simulated binary crossover
    Parent1 = [];
    Parent2 = [];
    for i=1:N/2
        Parent1 = [Parent1;population(MatingPool(i)).rnvec];
    end
    for i=N/2+1:N
        Parent2 = [Parent2;population(MatingPool(i)).rnvec];
    end
    % Parent1 = MatingPool(1:N/2,:); %First half of the mating population
    % Parent2 = MatingPool(N/2+1:end,:); 
    beta    = zeros(N/2*size(population(1).rnvec,1),D);
    miu     = rand(N/2*size(population(1).rnvec,1),D);
    beta(miu<=0.5) = (2*miu(miu<=0.5)).^(1/(DisC+1));
    beta(miu>0.5)  = (2-2*miu(miu>0.5)).^(-1/(DisC+1));
    beta = beta.*(-1).^randi([0,1],N/2*size(population(1).rnvec,1),D);
    beta(rand(N/2*size(population(1).rnvec,1),D)<0.5) = 1;
    beta(repmat(rand(N/2*size(population(1).rnvec,1),1)>ProC,1,D)) = 1;
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
%-----------------------------------------------------------------------------------------
% Polynomial mutation
    if rand<1 %Using the DTLZ mutation strategy
        MaxValue = repmat(Boundary(1,:),N*size(population(1).rnvec,1),1);
        MinValue = repmat(Boundary(2,:),N*size(population(1).rnvec,1),1);
        k    = rand(N*size(population(1).rnvec,1),D);
        miu  = rand(N*size(population(1).rnvec,1),D);
        Temp = k<=ProM & miu<0.5;
        Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(Offspring(Temp)-MinValue(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
        Temp = k<=ProM & miu>=0.5; 
        Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*(1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-Offspring(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));
    %-----------------------------------------------------------------------------------------
    % Set the offsprings which are infeasible to boundary values
        % Offspring(Offspring>MaxValue) = MaxValue(Offspring>MaxValue);
        % Offspring(Offspring<MinValue) = MinValue(Offspring<MinValue);
    
        newpop = population;
        for pos = 1:N
            % disp('DTLZ Mut') %DTLZ Mut is facing an issue as it is getting stuck in the adjust_constraint_turning_angle function
            % disp(pos)
            newpop(pos).rnvec = Offspring((pos-1)*size(population(1).rnvec,1)+1:pos*size(population(1).rnvec,1),:);
            for k = 1:size(newpop(pos).rnvec)
                newpop(pos).rnvec = Evolve.check_boundary(newpop(pos).rnvec,k,model);
            end
            newpop(pos).rnvec = sortrows(newpop(pos).rnvec,randi(2));
            newpop(pos).path = testBspline([newpop(pos).rnvec(:,1)';newpop(pos).rnvec(:,2)'],model.xmax)';            
            newpop(pos) = adjust_constraint_turning_angle(newpop(pos),model);            
            % newpop(pos).rnvec = Evolve.check_boundary(newpop(pos).rnvec,pos,model);
            newpop(pos) = evaluate(newpop(pos));
        end
    else %Using no mutation strategy, just crossover
        newpop = population;
        for pos = 1:N
            disp('No Mut')
            newpop(pos).rnvec = Offspring((pos-1)*size(population(1).rnvec,1)+1:pos*size(population(1).rnvec,1),:);
            for k = 1:size(newpop(pos).rnvec)
                newpop(pos).rnvec = Evolve.check_boundary(newpop(pos).rnvec,k,model);
            end
            newpop(pos).rnvec = sortrows(newpop(pos).rnvec,randi(2));
            newpop(pos).path = testBspline([newpop(pos).rnvec(:,1)';newpop(pos).rnvec(:,2)'],model.xmax)';
            % newpop(pos).rnvec = Evolve.check_boundary(newpop(pos).rnvec,pos,model);
            newpop(pos) = adjust_constraint_turning_angle(newpop(pos),model);
            newpop(pos) = evaluate(newpop(pos));
        end
    end


end
