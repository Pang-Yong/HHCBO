classdef HHCBO < ALGORITHM
% <multi> <real> <expensive>

% wmax --- 10 --- Number of generations before updating Kriging models
% ref  --- [1.1 1.1] --- reference point in normalized objective space
% ineObjIndex --- [1] --- index of inexpensive objectives
% ineConIndex --- [2] --- index of inexpensive constraints


%------------------------------- Reference --------------------------------
% Pang, Yong, et al. "Surrogate-assisted expensive constrained Bi-objective 
% optimization with highly heterogeneous evaluations." Swarm and Evolutionary 
% Computation 83 (2023): 101401.

% Author: Yong Pang
% pangy@mail.dlut.edu.cn



    methods
        function main(Algorithm,Problem)
            warning on backtrace
            warning off verbose
            %% Parameter setting
            [wmax,ref,ineObjIndex,ineConIndex] = Algorithm.ParameterSet(10,[1.1 1.1],[2],[1]);
            %% Initialization of NSGAII
            
            NI            = Problem.N;
            P             = UniformPoint(NI,Problem.D,'Latin');
            Population    = SOLUTION(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            A             = Population;
            

            
          %%  distinguish expensive and inexpensive
            
            eObjIndex=setdiff(1:Problem.M,ineObjIndex)   ;         
            CM=size(Population.cons,2)         ;  
            eConIndex=setdiff(1:CM,ineConIndex);
            
            if Problem.M>2
                error('Number of objective is larger than 2')
            elseif isempty(eObjIndex)
                error('No expensive objective')
            end

            thetaO         = 10.*ones(size(eObjIndex,2),Problem.D);
            thetaC         = 10.*ones(size(eConIndex,2),Problem.D);
            
            ModelO    = cell(1,size(eObjIndex,2));
            ModelC    = cell(1,size(eConIndex,2));         
            
            
            % main loop           
            while Algorithm.NotTerminated(A)
                Dec = Population.decs;
                Obj = Population.objs;
                Con = Population.cons;
                MSEO = zeros(size(Dec,1),size(eObjIndex,2));
                MSEC = zeros(size(Dec,1),size(eConIndex,2));
                
                
                %train surrogate model
                train_X = A.decs;
                YO = A.objs;
                YC = A.cons;
                train_YO = YO(:,eObjIndex);
                train_YC = YC(:,eConIndex);
                %size(train_X ,1)
                [~,distinct] = unique(round(train_X*1e6)/1e6,'rows');  
                train_X    = train_X(distinct,:);
                train_YO   = train_YO(distinct,:);
                train_YC   = train_YC(distinct,:);
                for i = 1:size(eObjIndex,2) % train objective surrogates
                    dmodel     = dacefit(train_X,train_YO(:,i),'regpoly0','corrgauss',thetaO(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    ModelO{i}   = dmodel;
                    thetaO(i,:) = dmodel.theta;
                end

                for i = 1:1:size(eConIndex,2) % train constraint surrogates
                    dmodel     = dacefit(train_X,train_YC(:,i),'regpoly0','corrgauss',thetaC(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    ModelC{i}   = dmodel;
                    thetaC(i,:) = dmodel.theta;
                end                
                
                
                % calculate real pareto front
                [RealFrontNo,~] = NDSort(YO,1);
                RealFirstObj = YO((RealFrontNo==1),:);
                
                
                
                [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
                w = 1;
                while w <= wmax
                    w = w + 1;
                    MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                    OffspringDec = OperatorGA(Dec(MatingPool,:));
                    N = size(OffspringDec,1);
                    OffspringObj = Problem.CalObj(OffspringDec);
                    OffspringCon = Problem.CalCon(OffspringDec);
                    OffspringMSEO = zeros(N,size(eObjIndex,2));
                    OffspringMSEC = zeros(N,size(eConIndex,2));

                    for i = 1:N                     
                        for j = 1:size(eObjIndex,2)
                            index=eObjIndex(j);
                            [OffspringObj(i,index),~,OffspringMSEO(i,j)] = predictor(OffspringDec(i,:),ModelO{j});
                        end
                        for j = 1:size(eConIndex,2)
                            index=eConIndex(j);
                            [OffspringCon(i,index),~,OffspringMSEC(i,j)] = predictor(OffspringDec(i,:),ModelC{j});
                        end
                    end                   
                    
                    all_Obj = [Obj;OffspringObj];
                    all_Con = [Con;OffspringCon];
                    all_MSEO = [MSEO;OffspringMSEO];
                    all_MSEC = [MSEC;OffspringMSEC];
                    all_Dec = [Dec;OffspringDec];
                    
                    eCon=all_Con(:,eConIndex);
                    ineCon=all_Con(:, ineConIndex);  
                    [Choose,FrontNo,CrowdDis] = EnvironmentalSelection2(all_Obj,eCon,ineCon,NI);
                    Dec = all_Dec(Choose,:);
                    Obj = all_Obj(Choose,:);
                    Con = all_Con(Choose,:);
                    MSEO = all_MSEO(Choose,:);
                    MSEC = all_MSEC(Choose,:);
                end
                
                
               %% Effienct EI sampling criterion 
                        
                EI=CalCEHVI(RealFirstObj,Obj,Con,eObjIndex,eConIndex,ineObjIndex,ineConIndex,MSEO,MSEC,ref);

                [~,sortIndex] = sort(EI,'descend');
                ChoseMax=max([round(size(EI,1)*unifrnd(0,0.3)),1]);
                FnewIndex=sortIndex(1:ChoseMax);   
                
                DA=RealFirstObj;
                dist_D=zeros(size(FnewIndex,1),size(DA,1));
                for i = 1:  size(FnewIndex,1 )
                    for j = 1:size(DA,1)
                        dist_D(i,j) = norm(Obj(FnewIndex(i,1),:)-DA(j,:),2);
                    end
                end
                
                % Diversity Indicator
                DI = min(dist_D,[],2);  
                [ ~,SnewIndex] =max(DI);
                newIndex=[FnewIndex(SnewIndex,1)];
                PnewDec = Dec(newIndex,:);
                
                
                PnewDec = unique(PnewDec,'rows');
                
                New = SOLUTION(PnewDec);

                A = [A,New];            
                [Population,~,~] = EnvironmentalSelection(A,Problem.N); 
            end      
        end
    end
end