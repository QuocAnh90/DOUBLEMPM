function [pore_wp] = Update_Pore_Pressure_new1(...
    wpCount,NODES,dt,active_elements,ElemsCount,nodeCount,mwpoints,CONNECT,nvelo_si,nvelo_wi,n_wp,n_i,pore_wp,dN,Kw)


de_wp = zeros(nodeCount,1);
I = ones(nodeCount,1);

for n = 1:nodeCount
 if norm(nvelo_si(n,:)) ~=0
nvelo_mixture(n,:) = (1-n_i(n))*nvelo_wi(n,:) - n_i(n)*nvelo_si(n,:); 
 else
 nvelo_mixture = nvelo_wi;
 end
end
 
L_mixture = cell(wpCount,1);
dWP   = cell(wpCount,1); 

% Calculate pore water pressure for liquid phase
for i=1:ElemsCount
    c = active_elements(i);                                                 % loop all the element
    mptw = mwpoints{c};                                                % all particles in the element
    
    for wp = 1:length(mptw)                                            % loop all particle in that element
        wpid = mptw(wp);                                               % wpid: global index of the particle
        L_mixture{wpid} = zeros(2,2);                                  % Strain rate tensor of liquid phase
    
for j=1:NODES(wpid)                                                    % loop all nodes interacting with particle wpid
          if dN{wpid}(j)==0
         continue
          end
            npid = CONNECT{wpid}(j);                                   % npid: global index of node
            L_mixture{wpid} = L_mixture{wpid} + (dN{wpid}(:,j)*nvelo_mixture(npid,:));          % Interpolate the strain rate from nodes to particles
end
        
        dWP{wpid} = (L_mixture{wpid} + L_mixture{wpid}')/2;      
        de_wp(wpid)  = dt/n_wp(wpid) * trace(dWP{wpid});
          
        pore_wp(wpid) = pore_wp(wpid) + Kw * de_wp(wpid);
%         pore_wp(wpid,1) = pore_wp(wpid,1) + Kw * de_wp(wpid);
%         pore_wp(wpid,2) = pore_wp(wpid,2) + Kw * de_wp(wpid);
    end      
end