function [Score] = calMetirc(MetricIndex,PopObj,problemIndex)

if MetricIndex == 1
    if problemIndex == 1
        a = [4000,4000];
    elseif problemIndex == 2
        a = [50,10];
    else
        a = [4000,4000];
        %disp('hello')
    end
    Score = HV(PopObj,a);
else
    Score = PD(PopObj);
end

end

