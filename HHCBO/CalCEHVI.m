function CEI = CalCEHVI(RealFirstObj,Obj,Con,eObjIndex,eConIndex,ineObjIndex,ineConIndex,MSEO,MSEC,ref)

%EHVI exact calculation 


      NRealFront = size(RealFirstObj,1);
       
      PopObj=[RealFirstObj;Obj];
      [NPop,M] = size(Obj);
      Zmin=min(PopObj,[],1);
      Zmax=max(PopObj,[],1);
      
      %% Normalization
      a=Zmax-Zmin;
      % Normalization      
 
      
      delta = 1;
      Obj = Obj - repmat(Zmin,size(Obj,1),1);
      RealFirstObj = RealFirstObj - repmat(Zmin,size(RealFirstObj,1),1);
      Obj = Obj./repmat(a,NPop,1);
      RealFirstObj  = RealFirstObj ./repmat(a,NRealFront,1);
      
      MSEO=MSEO./repmat((a(eObjIndex).^2),size(MSEO,1),1);
      sigmaO=MSEO.^0.5;
      sigmaC=MSEC.^0.5;
      
      
      
      if isempty(ineObjIndex)
          % two objectives are expensive
          EI=zeros(NPop,1);
          CEI=zeros(NPop,1);
          RealFirstObj = sortrows(RealFirstObj,1);
          aguPF=[-inf,ref(2); RealFirstObj; ref(2),-inf];
          a=aguPF(:,1);
          b=aguPF(:,2);          
  
          
           for i =1: NPop  
              if any(sigmaO(i,:)==0) | any(Con(i,ineConIndex)>0)
                  EI(i)=0;
              else
                  for k=1:(size(RealFirstObj,1))
                      for j=1:k-1
                            EI(i)=EI(i)+(b(j)-b(j+1))*(Pi(b(k),Obj(i,2),sigmaO(i,2))-Pi(b(k+1),Obj(i,2),sigmaO(i,2)))*Ei(a(j+1),Obj(i,1),sigmaO(i,1));
                      end
                      EI(i)=EI(i)+(Ei(b(k),Obj(i,2),sigmaO(i,2))-(b(k)-b(k+1))*Pi(b(k+1),Obj(i,2),sigmaO(i,2))-Ei(b(k+1),Obj(i,2),sigmaO(i,2)) )*Ei(a(k+1),Obj(i,1),sigmaO(i,1));  
                  end
                  k=size(RealFirstObj,1)+1;
                  for j=1:k-1
                      EI(i)=EI(i)+(b(j)-b(j+1))*Pi(b(k),Obj(i,2),sigmaO(i,2))*Ei(a(j+1),Obj(i,1),sigmaO(i,1));
                  end
                  EI(i)=EI(i)+Ei(b(k),Obj(i,2),sigmaO(i,2))*Ei(a(k+1),Obj(i,1),sigmaO(i,1));
                  
                  
                  %%POF
                  PF = 1;
                  if size(eConIndex,2) ~= 0
                      for j =  size(eConIndex,2)
                          index=eConIndex(j);
                          pf=normcdf(0 ,Con(i,index),sigmaC(i,j));
                          modpf=min(delta,pf);
                          PF = PF*modpf;
                      end
                  end
                   CEI(i)=EI(i)*PF;
              
              end                          
           end
           
      else          
          % one objectives are expensive    
          EI=zeros(NPop,1);
          CEI=zeros(NPop,1);
          RealFirstObj = [RealFirstObj(:,eObjIndex),RealFirstObj(:,ineObjIndex)];
          RealFirstObj = sortrows(RealFirstObj,1);
          aguPF=[-inf,ref(2); RealFirstObj; ref(2),-inf];
          a=aguPF(:,1);
          b=aguPF(:,2);
          for i =1: NPop                 
              if any(sigmaO(i,:)==0) | any(Con(i,ineConIndex)>0)
                  EI(i)=0;
              else
                  List = find(b>Obj(i,ineObjIndex));
                  I = List(end);
                  for j=1:I-1
                      L=a(j+1)-Obj(i,eObjIndex);
                      EI(i)=EI(i)+((b(j)-b(j+1))*(L*normcdf(L/sigmaO(i,1))+sigmaO(i,1)*normpdf(L/sigmaO(i,1))));
                  end    
                  L=a(I+1)-Obj(i,eObjIndex);
                  EI(i)=EI(i)+((b(I)-Obj(i,ineObjIndex))*(L*normcdf(L/sigmaO(i,1))+sigmaO(i,1)*normpdf(L/sigmaO(i,1))));

              
                  %%POF
                  PF = 1;
                  if size(eConIndex,2) ~= 0
                      for j =  size(eConIndex,2)
                          index=eConIndex(j);
                          pf=normcdf(0 ,Con(i,index),sigmaC(i,j));
                          modpf=min(delta,pf);
                          PF = PF*modpf;
                      end
                  end
                  CEI(i)=EI(i)*PF;
                     
              end                 
          end        
      end       
end      
      


function pi=Pi(a,mu,sigma)
    pi=normpdf(a,mu,sigma);
end


function ei=Ei(a,mu,sigma)
    ei=(a-mu)*normcdf((a-mu)/sigma)+sigma*normpdf((a-mu)/sigma);
end