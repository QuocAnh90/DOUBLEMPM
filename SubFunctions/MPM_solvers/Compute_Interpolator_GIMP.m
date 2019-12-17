function [N,dN,CONNECT,pElems,mpoints,NODES] = Compute_Interpolator_GIMP(pCount,cellCount,x_p,le,NN,LOC,lp)

pElems         = zeros(pCount,1);                             % index of elements where stores particles
CONNECT_TEMP    = zeros(pCount,16);                            % node 1=leftdown 2=righdown 3=rightup 4= leftup
NODES           = 16 * ones(pCount,1);                         % Number of interation nodes for each particle
N_local         = zeros(pCount,16);                            % Value of shape function
dN_local        = zeros(pCount,32);                            % Value of gradient of shape function

N = cell(pCount,1);
dN = cell(pCount,1);
CONNECT = cell(pCount,1);

 for p = 1:pCount
%  pElems(p) = ceil(x_p(p,1)/le(1))+(NN(1)-1)*(fix(x_p(p,2)/le(2)));   % compute vector store index elements 
 pElems(p) = floor(x_p(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_p(p,2)/le(2)));
 
 CONNECT_TEMP(p,6) = pElems(p) + ceil(pElems(p)/(NN(1)-1)) - 1;
 CONNECT_TEMP(p,7) = CONNECT_TEMP(p,6) + 1;
 CONNECT_TEMP(p,11)= CONNECT_TEMP(p,7) + NN(1);
 CONNECT_TEMP(p,10)= CONNECT_TEMP(p,6) + NN(1);
 CONNECT_TEMP(p,1) = CONNECT_TEMP(p,6) - NN(1) - 1;
 CONNECT_TEMP(p,2) = CONNECT_TEMP(p,1) + 1;
 CONNECT_TEMP(p,3) = CONNECT_TEMP(p,2) + 1;
 CONNECT_TEMP(p,4) = CONNECT_TEMP(p,3) + 1;
 CONNECT_TEMP(p,5) = CONNECT_TEMP(p,6) - 1;
 CONNECT_TEMP(p,8) = CONNECT_TEMP(p,7) + 1;
 CONNECT_TEMP(p,9) = CONNECT_TEMP(p,5) + NN(1);
 CONNECT_TEMP(p,12)= CONNECT_TEMP(p,11) + 1;
 CONNECT_TEMP(p,13)= CONNECT_TEMP(p,9) + NN(1);
 CONNECT_TEMP(p,14)= CONNECT_TEMP(p,13) + 1;
 CONNECT_TEMP(p,15)= CONNECT_TEMP(p,14) + 1;
 CONNECT_TEMP(p,16)= CONNECT_TEMP(p,15) + 1;
 
  for i=1:16
 CONNECT{p} = [CONNECT{p} CONNECT_TEMP(p,i)];
  end
 
for i = 1:NODES(p)
     % Compute the shape functions and gradient of the shape functions
    [N_local(p,i),dN_local(p,i),dN_local(p,i+8)]=GIMPshape(x_p(p,1:2),LOC(CONNECT_TEMP(p,i),:),le(1),le(2),lp{p});
    N{p}(i) = N_local(p,i);    
    dN{p}(1,i) = dN_local(p,i);
    dN{p}(2,i) = dN_local(p,i+8);
 end
 end
 
 % Compute mspoints: index of particles in each element (active element)
 for c =1:cellCount
     id_p = find(pElems==c);
     mpoints{c}=id_p;
 end
 