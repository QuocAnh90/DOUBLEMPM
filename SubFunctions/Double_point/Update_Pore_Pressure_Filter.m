function [pore_wp] = Update_Pore_Pressure_Filter(wpCount,NODES,dt,cellCount,mwpoints,CONNECT,nvelo_si,nvelo_wi,n_wp,pore_wp,dN,Kw)

de_wp  = zeros(wpCount,1);
de_swp = zeros(wpCount,1);
de_wwp = zeros(wpCount,1);
L_wp   = cell(wpCount,1);
L_swp  = cell(wpCount,1);
dEWP   = cell(wpCount,1); 
dESWP  = cell(wpCount,1);

% Calculate pore water pressure for liquid phase
for c=1:cellCount                                                      % loop all the element
    mptw = mwpoints{c};                                                % all particles in the element

    if isempty(mptw) ==0
    dEWP_F = zeros(2,2*length(mptw));
    dS = zeros(2*length(mptw),4);
    
    for wp = 1:length(mptw)                                            % loop all particle in that element
        wpid = mptw(wp);                                               % wpid: global index of the particle
        L_wp{wpid} = zeros(2,2);                                             % Strain rate tensor of liquid phase
        L_swp{wpid} = zeros(2,2);                                            % Strain rate tensor of solid phase
    
    for j=1:NODES(wpid)                                                    % loop all nodes interacting with particle wpid
          if dN{wpid}(j)==0
         continue
          end
            npid = CONNECT{wpid}(j);                                   % npid: global index of node
            L_wp{wpid} = L_wp{wpid} + (nvelo_wi(npid,:)'*dN{wpid}(:,j)');          % Interpolate the strain rate from nodes to particles
            L_swp{wpid} = L_swp{wpid} + (nvelo_si(npid,:)'*dN{wpid}(:,j)');
    end
        
        dEWP{wpid} = (L_wp{wpid} + L_wp{wpid}')/2; 
        dESWP{wpid} = (L_swp{wpid} + L_swp{wpid}')/2;
             
        de_wp(wpid)  = dt/n_wp(wpid) * (1-n_wp(wpid)) * trace(dEWP{wpid});
        de_swp(wpid) = dt/n_wp(wpid) * (-n_wp(wpid)) * trace(dESWP{wpid});      
        de_wwp(wpid) = de_wp(wpid) + de_swp(wpid);
     
        dEWP_F(:,wp*2-1:wp*2) = [de_wwp(wpid) 0;0 de_wwp(wpid)];        
        dS(2*wp-1:2*wp,:) = dN{wpid};
    end      

    % QR method    
        [Q1,R]=qr(dS,0);        
        r = rank(R);
%         if r~=3
%             break
%         end
        dEWP_F = (dEWP_F*Q1(:,1:r))*Q1(:,1:r)';
    end

    % Update pore water pressure
    for wp = 1:length(mptw)
        wpid = mptw(wp); 

        PORE = dEWP_F(:,wp*2-1:wp*2);
        de_wwp(wpid) = trace(PORE);
         pore_wp(wpid)  = pore_wp(wpid) + Kw * de_wwp(wpid);
    end
end
