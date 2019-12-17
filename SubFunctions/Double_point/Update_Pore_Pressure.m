function [pore_wp] = Update_Pore_Pressure(wpCount,NODES,dt,cellCount,mwpoints,CONNECT,nvelo_si,nvelo_wi,n_wp,pore_wp,dN,Kw)

de_wp  = zeros(wpCount,1);
de_swp = zeros(wpCount,1);
de_nwp = zeros(wpCount,1);
de_wwp = zeros(wpCount,1);
L_wp   = cell(wpCount,1);
L_swp  = cell(wpCount,1);
dEWP   = cell(wpCount,1); 
dESWP  = cell(wpCount,1);

% Calculate pore water pressure for liquid phase
for c=1:cellCount                                                      % loop all the element
    mptw = mwpoints{c};                                                % all particles in the element
    
    for wp = 1:length(mptw)                                            % loop all particle in that element
        wpid = mptw(wp);                                               % wpid: global index of the particle
        L_wp{wpid} = zeros(2,2);                                             % Strain rate tensor of liquid phase
        L_swp{wpid} = zeros(2,2);                                            % Strain rate tensor of solid phase
    
for j=1:NODES(wpid)                                                    % loop all nodes interacting with particle wpid
          if dN{wpid}(j)==0
         continue
          end
            npid = CONNECT{wpid}(j);                                   % npid: global index of node
            L_wp{wpid} = L_wp{wpid} + (dN{wpid}(:,j)*nvelo_wi(npid,:));          % Interpolate the strain rate from nodes to particles
            L_swp{wpid} = L_swp{wpid} + (dN{wpid}(:,j)*nvelo_si(npid,:));
end
        
        dEWP{wpid} = (L_wp{wpid} + L_wp{wpid}')/2; 
        dESWP{wpid} = (L_swp{wpid} + L_swp{wpid}')/2;
        
        de_wp(wpid)  = dt/n_wp(wpid) * (1-n_wp(wpid)) * trace(dEWP{wpid});
        de_swp(wpid) = dt/n_wp(wpid) * (-n_wp(wpid)) * trace(dESWP{wpid});
%         de_nwp(wpid) = dt/n_wp(wpid) * ((nvelo_wi(npid,:)-nvelo_si(npid,:)) * dn_wp);
        
%         de_wp(wpid)  = dt/n_wp(wpid) * (n_wp(wpid)) * trace(dEWP{wpid});
%         de_swp(wpid) = dt/n_wp(wpid) * (1-n_wp(wpid)) * trace(dESWP{wpid});
        
        de_wwp(wpid) = de_wp(wpid) + de_swp(wpid) + de_nwp(wpid);
        
        pore_wp(wpid) = pore_wp(wpid) + Kw * de_wwp(wpid);
        
    end      
end