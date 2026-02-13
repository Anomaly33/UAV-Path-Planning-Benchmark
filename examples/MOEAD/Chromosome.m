 classdef Chromosome
    
    properties
        rnvec;
        path;
        objs;
        front;
        vel;
        CD;
        rank;
        cons;
        dominationcount=0;
        dominatedset=[];
        dominatedsetlength=0;
        highBound;
        Qtable;
        currentState;
        reward;
    end
    
    methods

        function object=Chromosome(model)
            dim = model.n;
            bound = model.xmax;
            object.highBound(1) = model.xmin;
            for i = 2 : dim-1
                object.highBound(i) = i*bound/model.n; 
            end

            object.rnvec(1,:) = model.start; 
            object.rnvec(model.n,:) = model.end;
        end

        function object=initialize(object,model)
            if model.ymax == 200
                varRange = [-20,20];
            else
                varRange = [-5,5];
            end
            bound = model.xmax;
            i = 2;
            while i < model.n
                object.rnvec(i,2) = i*bound/model.n + rand(1)*(varRange(2)-varRange(1))+varRange(1);
                object.rnvec(i,1) = i*bound/model.n + rand(1)*(varRange(2)-varRange(1))+varRange(1);
                if ~(object.rnvec(i,1) < object.rnvec(i-1,1) || object.rnvec(i,2) < object.rnvec(i-1,2))
                    i = i + 1;
                end
            end
            
            a = object.rnvec(:,1);
            a(a < model.ymin) = model.xmin; 
            a(a > model.ymax) = model.ymax;
            object.rnvec(:,1) = a;
            a = object.rnvec(:,2);
            a(a < model.ymin) = model.xmin; 
            a(a > model.ymax) = model.ymax;
            object.rnvec(:,2) = a;
            object.path = testBspline([object.rnvec(:,1)';object.rnvec(:,2)'],model.xmax)'; 
            
            object = adjust_constraint_turning_angle(object,model);

        end
        
        function object = custom_cv(object,model)
            const_viol = 0;
            for k = 2 : size(object.path,1)-1 
                if k>3
                    L1 = sqrt((object.path(k,1)-object.path(k-1,1))^2+(object.path(k,2)-object.path(k-1,2))^2);
                    L2 = sqrt((object.path(k-1,1)-object.path(k-2,1))^2+(object.path(k-1,2)-object.path(k-2,2))^2);
                    L3 = sqrt((object.path(k,1)-object.path(k-2,1))^2+(object.path(k,2)-object.path(k-2,2))^2);
                    alpha = acosd((L1^2+L2^2-L3^2)/(2*L1*L2));
                    if alpha < 75
                        const_viol = const_viol+abs(alpha-75);
                    end
                
                    if object.path(k,1) < model.xmin
                        object.path(k,1) = model.xmin;
                    end
                    if object.path(k,1) > model.xmax
                        object.path(k,1) = model.xmax;
                    end
                    if object.path(k,2) < model.ymin
                        object.path(k,2) = model.ymin;
                    end
                    if object.path(k,2) > model.ymax
                        object.path(k,2) = model.ymax;
                    end
                end
                L4 = sqrt((object.path(k,1)-object.path(k-1,1))^2+(object.path(k,2)-object.path(k-1,2))^2); 
                beta = atand(abs(object.path(k,3)-object.path(k-1,3))/L4); 
                if beta > 60
                    const_viol = const_viol + abs(beta-60);
                end
            end
            object.cons = const_viol;
        end
    
        function object = get_cv(object,model)            
            cv = 0;
            object.path(1,3) = model.start(3);%+1500;
            for i = 2 : size(object.path,1)-1                 
                if i>3
                    L1 = sqrt((object.path(i,1)-object.path(i-1,1))^2+(object.path(i,2)-object.path(i-1,2))^2);
                    L2 = sqrt((object.path(i-1,1)-object.path(i-2,1))^2+(object.path(i-1,2)-object.path(i-2,2))^2);
                    L3 = sqrt((object.path(i,1)-object.path(i-2,1))^2+(object.path(i,2)-object.path(i-2,2))^2);
                    alpha = acosd((L1^2+L2^2-L3^2)/(2*L1*L2));
                    if alpha < 75
                        cv = cv+abs(alpha-75);
                    end
                
                    if object.path(i,1) < model.xmin
                        object.path(i,1) = model.xmin;
                    end
                    if object.path(i,1) > model.xmax
                        object.path(i,1) = model.xmax;
                    end
                    if object.path(i,2) < model.ymin
                        object.path(i,2) = model.ymin;
                    end
                    if object.path(i,2) > model.ymax
                        object.path(i,2) = model.ymax;
                    end
                    if object.path(i,3)<model.zmin
                        object.path(i,3) = model.zmin;
                    end
                    if object.path(i,3)>model.zmax
                        object.path(i,3) = model.zmax;
                    end
                end
                L4 = sqrt((object.path(i,1)-object.path(i-1,1))^2+(object.path(i,2)-object.path(i-1,2))^2); 
                beta = atand(abs(object.path(i,3)-object.path(i-1,3))/L4); 
                if beta > 60
                    cv = cv + abs(beta-60);
                end

                if model.xmax == 20
                    v1 = [model.H(floor(object.path(i,2))*10,floor(object.path(i,1))*10);   
                    model.H(floor(object.path(i,2))*10,ceil(object.path(i,1))*10);
                    model.H(ceil(object.path(i,2))*10,floor(object.path(i,1))*10);
                    model.H(ceil(object.path(i,2))*10,ceil(object.path(i,1))*10)];
                else
                    v1 = [model.H(floor(object.path(i,2)),floor(object.path(i,1)));
                    model.H(floor(object.path(i,2)),ceil(object.path(i,1)));
                    model.H(ceil(object.path(i,2)),floor(object.path(i,1)));
                    model.H(ceil(object.path(i,2)),ceil(object.path(i,1)))];
                end
                object.path(i,3) = max(v1) + model.safeH;
            end
            object.path(end,3) = model.end(3);
            object.cons = cv;
        end
                

        function [object] = adjust_constraint_turning_angle(object,model)
            
            object.path(1,3) = model.H(1,1);
            object.path(1,1) = model.start(1);
            object.path(1,2) = model.start(2);
            if model.ymax == 200
                varRange = [-2,10];
            else
                varRange = [-5,5];
            end
            for i = 2 : size(object.path,1)-1
                object.path = Evolve.check_boundary(object.path,i,model);
                while i > 3 && check_constraint_horizontal_turning_angle(object,i) ~= 0 
                    if rand(1) < 0.5
                        object.path(i,2) = object.path(i-1,2) + rand*(varRange(2)-varRange(1))+varRange(1);
                    else
                        object.path(i,1) = object.path(i-1,1) + rand*(varRange(2)-varRange(1))+varRange(1);                                                                        
                    end
                    object.path = Evolve.check_boundary(object.path,i,model); 
                end
                if model.xmax == 20
                    v1 = [model.H(floor(object.path(i,2))*10,floor(object.path(i,1))*10);   
                    model.H(floor(object.path(i,2))*10,ceil(object.path(i,1))*10);
                    model.H(ceil(object.path(i,2))*10,floor(object.path(i,1))*10);
                    model.H(ceil(object.path(i,2))*10,ceil(object.path(i,1))*10)];
                else
                    v1 = [model.H(floor(object.path(i,2)),floor(object.path(i,1)));
                    model.H(floor(object.path(i,2)),min(model.xmax,ceil(object.path(i,1))));
                    model.H(min(model.ymax,ceil(object.path(i,2))),floor(object.path(i,1)));
                    model.H(min(model.ymax,ceil(object.path(i,2))),min(model.xmax,ceil(object.path(i,1))))];
                end

                    object.path(i,3) = max(v1) + model.safeH; 
                while check_constraint_vertical_turning_angle(object,i) ~= 0 
                    j = i; 
                    while j > 1
                        if object.path(j,3) < object.path(j-1,3) 
                            object.path(j,3) = object.path(j,3) + rand*(object.path(j-1,3)-object.path(j,3));
                            if check_constraint_vertical_turning_angle(object,j) == 0
                                break;
                            end
                        else
                            object.path(j-1,3) = object.path(j-1,3) + rand*(object.path(j,3)-object.path(j-1,3));
                            if check_constraint_vertical_turning_angle(object,j) == 0
                                if j > 2 && check_constraint_vertical_turning_angle(object,j-1) == 0
                                    break;
                                else
                                    j = j - 1;
                                end
                            end
                                
                        end
                    end                    
                end
            end
            object.path(end,3) = model.end(3);
            object.path(end,1) = model.end(1);
            object.path(end,2) = model.end(2);
            object.path(1,3) = model.H(1,1);
            object.path(1,1) = model.start(1);
            object.path(1,2) = model.start(2);
        end
        

        function [flag] = check_constraint_horizontal_turning_angle(object,i)
            flag = 0;
            L1 = sqrt((object.path(i,1)-object.path(i-1,1))^2+(object.path(i,2)-object.path(i-1,2))^2);
            L2 = sqrt((object.path(i-1,1)-object.path(i-2,1))^2+(object.path(i-1,2)-object.path(i-2,2))^2);
            L3 = sqrt((object.path(i,1)-object.path(i-2,1))^2+(object.path(i,2)-object.path(i-2,2))^2);
            alpha = acosd((L1^2+L2^2-L3^2)/(2*L1*L2));
            if alpha < 75
                flag = 1;
            end
        end

        function [flag] = check_constraint_vertical_turning_angle(object,i)
            flag = 0;
            L1 = sqrt((object.path(i,1)-object.path(i-1,1))^2+(object.path(i,2)-object.path(i-1,2))^2); 
            beta = atand(abs(object.path(i,3)-object.path(i-1,3))/L1); 
            if beta > 60

                flag = 1;
            end
        end
        

        function object=evaluate(object)
            
            % Input solution
            x_all = object.path(:,1);
            y_all = object.path(:,2);
            z_all = object.path(:,3);

            N = size(x_all,1); % Full path length
            %============================================
            % J1 - path length
            J1 = 0;
            
            for i = 1:N-1
                diff = [x_all(i+1) - x_all(i);y_all(i+1) - y_all(i);z_all(i+1) - z_all(i)];% Path Length
                J1 = J1 + norm(diff);
            end

           
            % J2 - Height
            J2 = sum(z_all)/ length(z_all);% Height of path segments

            %==============================================

            % Evaluation Function
            object.objs = [J1, J2];
        end
        
        
    end    
    
end

