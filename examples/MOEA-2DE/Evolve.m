classdef Evolve    
    methods(Static)

        function [Offspring] = dimExplore(Parent,dim,model,strongDim,F)
            Site = rand(1,dim-2) < F;
            Site = [0,Site,0];
            Site(strongDim==1)=1;
            Offspring = Parent;
            for i = 2 : dim-1
                if Site(i) == 0
                    continue;
                end
                eachBound = Offspring.highBound(i);
                llowerBound = model.xmin;
                if i == 2
                    lowerBound = Offspring.highBound(i);
                else
                    lowerBound = Offspring.highBound(i-1);
                end
                if eachBound>model.xmax
                    eachBound = model.xmax;
                end
                if rand(1)<0.5
                    Offspring.rnvec(i,1) = rand(1)*(lowerBound - llowerBound) + llowerBound;
                    Offspring.rnvec(i,2) = rand(1)*(eachBound - lowerBound) + lowerBound;
                else
                    Offspring.rnvec(i,2) = rand(1)*(lowerBound - llowerBound) + llowerBound;
                    Offspring.rnvec(i,1) = rand(1)*(eachBound - lowerBound) + lowerBound;
                end
            end

        end

        
        
        function [Offspring1,Offspring2] = binary_crossover(Parent1,Parent2,D,proC,i)
            Offspring1 = Parent1;
            Offspring2 = Parent2;
            % One point crossover
            N = 1;
            k = repmat(1:D-1,N,1) > repmat(randi(D-1,N,1),1,D-1);
            k(repmat(rand(N,1)>proC,1,D-1)) = false;
            temp = Offspring1.rnvec(2:D,i);
            temp1 = Parent2.rnvec(2:D,i);
            temp(k) = temp1(k);
            Offspring1.rnvec(2:D,i) = temp;
            temp = Offspring2.rnvec(2:D,i);
            temp1 = Parent1.rnvec(2:D,i);
            temp(k) = temp1(k);
            Offspring2.rnvec(2:D,i) = temp;
        end
    
        function [c] = mutation(p,dim,model)
            if model.ymax == 200
                varRange = [-10,10];
            else
                varRange = [-5,5];
            end
            pos = randi([2,dim-1],1,1);
            c = p;
            if rand < 0.5
                c.rnvec(pos,1) = pos*model.xmax/model.n + rand*(varRange(2)-varRange(1))+varRange(1);
            else
                c.rnvec(pos,2) = pos*model.ymax/model.n + rand*(varRange(2)-varRange(1))+varRange(1);
            end
            c.rnvec = Evolve.check_boundary(c.rnvec,pos,model);
        end

        function [p] = check_boundary(p,i,model)
            if p(i,1) < model.xmin
                p(i,1) = model.xmin;
            end
            if p(i,1) > model.xmax
                p(i,1) = model.xmax;
            end
            if p(i,2) < model.ymin
                p(i,2) = model.ymin;
            end
            if p(i,2) > model.ymax
                p(i,2) = model.ymax;
            end
        end

    end
end