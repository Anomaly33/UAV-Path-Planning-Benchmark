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
% figure
% %plot3(cpts(1,:),cpts(2,:),cpts(3,:),'xb-')
% plot(cpts(1,:),cpts(2,:),'xb-')
% hold all
% %plot3(q(1,:), q(2,:), q(3,:))
% plot3(q(1,:), q(2,:))
% xlabel('X')
% ylabel('Y')
% %zlabel('Z')
% hold off

end