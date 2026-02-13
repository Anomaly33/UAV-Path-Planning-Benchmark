function q = testBspline(cpts,n)

tpts = [0 5];
if n == 20
    m = 0.1;
else
%     m = 0.01;%500
    m = 0.025;%200
end
tvec = 0:m:5;
[q] = bsplinepolytraj(cpts,tpts,tvec);

end
