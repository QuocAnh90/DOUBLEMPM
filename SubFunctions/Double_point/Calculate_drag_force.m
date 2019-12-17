function [ndrag,niforce_wi] = Calculate_drag_force(...
    NODES,nodeCount,CONNECT,N,wpCount,m_wp,k_wp,nmass_wi,nmass_si,nmomentum_wi,nmomentum_si,alpha,niforce_cell)

drag                 = zeros(nodeCount,1); 
ndrag                = zeros(nodeCount,2); 
niforce_wi           = zeros(nodeCount,2);    
nvelo_wi             = zeros(nodeCount,2);
nvelo_si             = zeros(nodeCount,2);

%% Calculate the drag force
 for p=1:wpCount
 for j=1:NODES(p)
     npid                           = CONNECT{p}(j);
     
          if N{p}(j)==0
         continue
          end
     
% Drag forces
 drag(npid)          = drag(npid) + m_wp(p)*9.81*N{p}(j)/k_wp(p);
 end 
 end
 
 for n = 1:nodeCount
    if nmass_wi(n) == 0
        continue
    end
    nvelo_wi(n,:) = nmomentum_wi(n,:)/nmass_wi(n);
 end
 
  for n = 1:nodeCount
    if nmass_si(n) == 0
        continue
    end
    nvelo_si(n,:) = nmomentum_si(n,:)/nmass_si(n);
 end
 
for n = 1:nodeCount
    ndrag(n,:) = drag(n) * (nvelo_wi(n,:) - nvelo_si(n,:));
end

%% Calculate the internal force
 
for n=1:nodeCount
    niforce_wi(n,:)          = alpha(n) * niforce_cell(n,:);
end