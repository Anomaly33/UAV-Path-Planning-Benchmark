classdef BSpline < handle
    
    properties
        CVs
        nCV
        umin
        umax
        knots
        order
    end
    
    methods
        function idx = span_index(this,u)
            if u == this.umax
                idx = lower_bound(this.knots,u)-1;
            else
                idx = upper_bound(this.knots,u)-1;
            end
        end

        function idx = lower_bound(vector,x)
            count = length(vector);
            first = 1;
            while(count > 0)
                step = floor(count/2);
                it = first+step;
                if vector(it) < x
                    first = it+1;
                    count = count-(step+1);
                else
                    count = step;
                end
            end
            idx = first;
        end

        function idx = upper_bound(vector,x)
            count = length(vector);
            first = 1;
            while(count > 0)
                step = floor(count/2);
                it = first+step;
                if ~(vector(it) > x)
                    first = it+1;
                    count = count-(step+1);
                else
                    count = step;
                end
            end
            idx = first;
        end
    end
    
    methods
        function obj = BSpline(CVs,order,umin,umax)
            if nargin < 4
                umin = 0;
                umax = 1;
            end   
            nCV = size(CVs,1);
            if nCV < order
                obj = BSpline.empty;
                return
            end
            obj.CVs = CVs;
            obj.nCV = nCV;
            obj.order = order;
            knots = zeros(obj.nCV+obj.order,1);
            % left part has order knots
            for i = 1:order
                knots(i) = umin;
            end
            % middle part has nCV-order knots
            n_pnt = nCV-order;
            stride = (umax-umin)/(n_pnt+1);
            for i = 1:n_pnt
                knots(order+i) = umin+stride*i;
            end
            % right part has order knots
            for i = 1:order
                knots(order+n_pnt+i) = umax;
            end
            obj.umin = umin;
            obj.umax = umax;
            obj.knots = knots;
        end
        
        function this = set_knots(this,knots)
            if length(knots) ~= length(this.knots)
                return
            end
            if all(diff(knots) < 0)
                return
            end
            l_umin = min(knots);
            if ~all(knots(1:this.order) == l_umin)
                return
            end
            l_umax = max(knots);
            if ~all(knots(end-this.order+1:end) == l_umax)
                return
            end
            if ~all((knots(this.order+1:end-this.order) < l_umax) .* ...
                    (knots(this.order+1:end-this.order) > l_umin))
                return
            end
            this.umin = l_umin;
            this.umax = l_umax;
            this.knots = knots;
        end
        
        function N_i = basis(this,i,u)
            u = max(u,this.umin);
            u = min(u,this.umax);
            n = this.order-1;
            span_idx = span_index(this,u);
            if span_idx < i || span_idx > i+n
                N_i = 0;
                return
            end
            N_i_ns1 = zeros(1,this.order);
            N_i_n = zeros(1,this.order);
            N_i_ns1(span_idx-i+1) = 1;
            for nn = 1:n
                idx = 1;
                for ii = i:i+n-nn
                    N_a = N_i_ns1(idx);
                    N_b = N_i_ns1(idx+1);
                    if N_a ~= 0
                        c_a = coef(this,ii,nn,u);
                        N_a = c_a*N_a;
                    end
                    if N_b ~= 0
                        c_b = 1-coef(this,ii+1,nn,u);
                        N_b = c_b*N_b;
                    end
                    N_i_n(idx) = N_a+N_b;
                    idx = idx+1;
                end
                N_i_ns1 = N_i_n;
            end
            N_i = N_i_n(1);
            % support function
            function f_i_n = coef(obj,i,n,u)
                f_i_n = (u-obj.knots(i))/(obj.knots(i+n)-obj.knots(i));
            end
        end
        
        function [x,y] = evaluate(this,u)
            u = max(u,this.umin);
            u = min(u,this.umax);
            x = 0;
            y = 0;
            span_idx = span_index(this,u);
            n = this.order-1;
            for i = span_idx-n:span_idx
                N = this.basis(i,u);
                x = x+N*this.CVs(i,1);
                y = y+N*this.CVs(i,2);
            end
        end
        
        function [xs,ys] = evaluate_batch(this,us)
            xs = 0*us;
            ys = 0*us;
            for i = 1:length(us)
                [xs(i),ys(i)] = evaluate(this,us(i));
            end
        end

        
    end
    
end