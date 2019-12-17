function N_MLS = LinearNode_MLS(nodeCount,spCount,NODES,CONNECT_TEMP,CONNECT,x_sp,LOC,le,r,r_n)

% Funtion for mapping from particles to nodes using moving least square
% Least square matrix
A               = cell(nodeCount,1);

for i = 1:nodeCount
    A{i} = zeros(r_n(i),r_n(i));
end
% Funtion for mapping from particles to nodes using linear basis
N_local         = zeros(spCount,4);              % Value of shape function

for p = 1:spCount
for i = 1:NODES(p)
    npid = CONNECT_TEMP(p,i);
     % Compute the shape functions and gradient of the shape functions
    [N_local(p,i),~,~]=linearshape(x_sp(p,1:2),LOC(npid,:),le(1),le(2));
end
 end

% Compute A
for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];    
     elseif r==6
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
     end
         
     PPxy = Pxy * Pxy';
     
     for i = 1:NODES(p)
         npid = CONNECT{p}(i);
         A{npid} = A{npid} + N_local(p,i) * PPxy;  
     end
end

% % Check rank D
% for n=1:nodeCount
% %     if rank(A{n})<6
%     if rcond(A{n})<=9.999999e-10
%         r_n(n) = 3;
%     end
% end

% Compute shape function and gradients
 for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];   
     elseif r==6
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
     end
     
     for i = 1:NODES(p)
         npid = CONNECT{p}(i);
         if r==1
         pn = 1;
         elseif r==3
         pn = [1 ; LOC(npid,1) ; LOC(npid,2)];   
         elseif r==6
         pn = [1 ; LOC(npid,1) ; LOC(npid,2) ; LOC(npid,1)^2 ; LOC(npid,1)*LOC(npid,1) ; LOC(npid,1)^2];   
         end
         
         if r_n(npid)==6
         rn = A{npid} \ pn;
         N_MLS{p}(i) = rn' * N_local(p,i) * Pxy;
         
         elseif r_n(npid)==3
         rn = A{npid}(1:3,1:3) \ pn(1:3);
         N_MLS{p}(i) = rn' * N_local(p,i) * Pxy(1:3);       
         end        
     end
 end