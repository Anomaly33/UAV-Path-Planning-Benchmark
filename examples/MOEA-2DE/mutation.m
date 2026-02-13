function [c] = mutation(p,dim,model)

pos = randi(dim,1,1);
c = p;
rnd = randi([1,3]);
varRange = 4;
if rnd == 1
    c.rnvec(pos,1) = p.rnvec(pos - 1,1) + rand*varRange;
    if c.rnvec(pos,1) < 1
        c.rnvec(pos,1) = 1;
    end
    if c.rnvec(pos,1) > 200
        c.rnvec(pos,1) = 200;
    end 
else
    c.rnvec(pos,2) = p.rnvec(pos - 1,2) + rand*varRange;
    if c.rnvec(pos,2) < 1
        c.rnvec(pos,2) = 1;
    end
    if c.rnvec(pos,2) > 200
        c.rnvec(pos,2) = 200;
    end 
end

c.rnvec(pos,3) = model.H(c.rnvec(pos,1),c.rnvec(pos,2)) + model.safeH;

end