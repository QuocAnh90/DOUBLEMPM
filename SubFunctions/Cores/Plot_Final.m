function StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,s_sp,v_ssp,spCount,r1_sp,r2_sp)

%  x_corner1              = zeros(spCount,2);
%  x_corner2              = zeros(spCount,2);
%  x_corner3              = zeros(spCount,2);
%  x_corner4              = zeros(spCount,2);
%  
displacement = zeros(spCount,1);
velocity = zeros(spCount,1);
stress = zeros(spCount,1);
for sp=1:spCount
    displacement(sp) = sqrt(d_sp(sp,1)^2+d_sp(sp,2)^2);
    velocity(sp) = sqrt(v_ssp(sp,1)^2+v_ssp(sp,2)^2);
%     stress(sp) = norm(s_sp);
    stress(sp) = (s_sp(sp,2));
end
% 
% % max(velocity)
% 
%  for sp=1:spCount
%  x_corner1(sp,:) = x_sp(sp,:) - r1_sp(sp,:) - r2_sp(sp,:);      % Position of corner 1
%  x_corner2(sp,:) = x_sp(sp,:) + r1_sp(sp,:) - r2_sp(sp,:);
%  x_corner3(sp,:) = x_sp(sp,:) + r1_sp(sp,:) + r2_sp(sp,:);
%  x_corner4(sp,:) = x_sp(sp,:) - r1_sp(sp,:) + r2_sp(sp,:);
%  end

StressProfile1=figure;
% set(StressProfile1, 'visible','off');
sz = 5;
color = stress;
% color = b_sp(:,1);
% scatter(x_sp(:,1),x_sp(:,2),sz,color,'filled');
scatter(x_wp(:,1),x_wp(:,2),sz,color,'filled');

% hold on
% for sp=1:spCount
%     x_cor = [x_corner1(sp,1) x_corner2(sp,1) x_corner3(sp,1) x_corner4(sp,1) x_corner1(sp,1)];
%     y_cor = [x_corner1(sp,2) x_corner2(sp,2) x_corner3(sp,2) x_corner4(sp,2) x_corner1(sp,2)];
%     plot(x_cor,y_cor,'b')
%     hold on
% end

grid on
axis([0,max(LOC(:,1)),0,max(LOC(:,2))]);
set(gca,'xtick',[0:le(1):max(LOC(:,1))]);
set(gca,'ytick',[0:le(2):max(LOC(:,2))]);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]);
h=colorbar;
colormap(jet(256))
% set(h, 'ytick', [0:0.2:1.2]);
%     caxis([0 0.1]);


% for wp=1:wpCount
%     velocity1(wp) = sqrt(v_wp(wp,1)^2+v_wp(wp,2)^2);
%     stress1(sp) = pore_wp(wp);
% end