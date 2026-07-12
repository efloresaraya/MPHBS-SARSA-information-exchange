
% Hybrid EO and DE (collaborative & Sequential)
%Amal Moharam et al.,Neural Computing & Applications, 2024
%"Economically optimized heat exchanger design: a synergistic approach using
%differential evolution and equilibrium optimizer within an evolutionary 
%algorithm framework"
%Abstract: This study introduces the CP-EODE algorithm, a novel hybrid
% of the Equilibrium Optimizer (EO) ,and the Differential Evolution (DE)
% algorithm. It addresses EO's tendency towards premature convergence by
% enhancing its exploration capabilities by DE. 
%%
clear all;
clc;
%% Setting
Run_no=51;        % Number of Runs
nPop = 50 ;       % Population Size
func_num=29;   %select function number,In cec2017:func_num can be from 1 to 30 
ObjectiveFunction=str2func('cec17_func');        % Objective Function
  
Dim = 30;            % Number of Dimensions
DimSize = [1 Dim];   
LB = -100;        %  Lower Bound
UB = 100;         %  Upper Bound
MaxIter = 6000;   % Maximum Number of Iterations
                  %so FunctionEvaluations=6000*nPop=300,000=Dim*10,000
                  
sAveFit=zeros(Run_no,MaxIter);  %Average Fitness

%% EO,DE parameters
V=1;      
a1=2;
a2=1;
GP=0.5;
FD = 0.5 ;
CR=0.9;

%% 
for irun=1:Run_no
disp(irun)
%initialize
% 4 best individuals
Ceq1=zeros(1,Dim);   Ceq1_fit=inf; 
Ceq2=zeros(1,Dim);   Ceq2_fit=inf; 
Ceq3=zeros(1,Dim);   Ceq3_fit=inf; 
Ceq4=zeros(1,Dim);   Ceq4_fit=inf;

% Empty Individual Structure
empty_individual.Position = [];
empty_individual.Fitness = [];
empty_individual.Covergence_curve = [];
empty_individual.fit_old = [];
empty_individual.Position_old = [];

% Initialize Population Array
pop = repmat(empty_individual, nPop, 1);

%% ------------ Initialize Population randomly
for i = 1:nPop
    pop(i).Position = unifrnd(LB, UB, DimSize);
end
%%   
for it = 1:MaxIter
    %disp(it)
    sumFit=0; % variable used for calculate Average Fitness
%% ---------- Bound Constraint
        for i=1:nPop  
            Flag4ub=pop(i).Position>UB;
            Flag4lb=pop(i).Position<LB;
            pop(i).Position=(pop(i).Position).*(~(Flag4ub+Flag4lb))+UB.*Flag4ub+UB.*Flag4lb; 
           % disp (pop(i).Position)
        end
%% ---------- Evaluate fitness 
 %1- Convert structure to cell
        for i=1:nPop
            vals(i,:) = [pop(i).Position];
        end
        %disp(vals)
            Fit = feval(ObjectiveFunction,vals',func_num);
        for i = 1:nPop
            pop(i).Fitness = Fit(i);
        end
                
 %% --------------- Select 4 best individuals from the population
       for i=1:nPop  
            if pop(i).Fitness<Ceq1_fit 
                  Ceq1_fit=pop(i).Fitness;  Ceq1=pop(i).Position;
            elseif pop(i).Fitness>Ceq1_fit && pop(i).Fitness<Ceq2_fit  
                  Ceq2_fit=pop(i).Fitness;  Ceq2=pop(i).Position;
            elseif pop(i).Fitness>Ceq1_fit && pop(i).Fitness>Ceq2_fit && pop(i).Fitness<Ceq3_fit
                  Ceq3_fit=pop(i).Fitness;  Ceq3=pop(i).Position;
            elseif pop(i).Fitness>Ceq1_fit && pop(i).Fitness>Ceq2_fit && pop(i).Fitness>Ceq3_fit && pop(i).Fitness<Ceq4_fit
                  Ceq4_fit=pop(i).Fitness;  Ceq4=pop(i).Position;
            end
       end
%% ---------------- Memory saving   
      if it==1
          for i=1:nPop
         pop(i).fit_old=pop(i).Fitness;  pop(i).Position_old=pop(i).Position;
          end
      end
       for i=1:nPop
         if pop(i).fit_old<pop(i).Fitness;
             pop(i).Fitness=pop(i).fit_old; pop(i).Position=pop(i).Position_old;
         end
       end
      %save new minimum 
       for i=1:nPop
         pop(i).Position_old=pop(i).Position;  pop(i).fit_old=pop(i).Fitness;
       end  
      %Save average fitness
       for i=1:nPop
           sumFit=sumFit+pop(i).Fitness;   
       end
       sAveFit(irun,it)=sumFit/nPop;
%% EO setting
      Ceq_ave=(Ceq1+Ceq2+Ceq3+Ceq4)/4;                    % averaged individuals 
      C_pool=[Ceq1; Ceq2; Ceq3; Ceq4; Ceq_ave];           % Equilibrium pool
      t=(1-it/MaxIter)^(a2*it/MaxIter)   ;            
      %% ---------- DE
    
     if (rem(it,2)==0)  %for sequential opertion bt EO and DE 
          %disp('DE')
       for i = 1:nPop
         %% --------- Pick indexes for random difference vector
           rv1 = floor(rand()* nPop) + 1; 
           while (rv1==i)
             rv1 = floor(rand()* nPop) + 1;
           end
             rv2 = floor(rand()* nPop) + 1;
           while ((rv2==rv1)||(rv2==i))
             rv2 = floor(rand()* nPop) + 1;
           end
             rv3 = floor(rand()* nPop) + 1;
           while ((rv3==rv2)||(rv3==rv1)||(rv3==i))
             rv3 = floor(rand()* nPop) + 1;
           end                  
       %% ------- Create noisy random vector
       nrv=pop(rv3).Position+FD*(pop(rv1).Position-pop(rv2).Position);
       % boundary constraints 
        for idx=1:Dim
           if (nrv(idx) > UB)
               nrv(idx) = UB ;
           elseif (nrv(idx) < LB)
                   nrv(idx) =LB ;
           end
        end
      %% Create trial vector use binomial crossover and 1  difference vector 
           rndnum = rand(1,Dim) < CR ;   % all random numbers < CR are 1 &( 0 otherwise)
           invmask = rndnum < 0.5   ;     % inverse mask to rndnum
           Trialv = invmask.*pop(i).Position + rndnum.*nrv   ;   % crossover
           pop(i).Position=Trialv;
       end
      %%  --------------- EO
     else  %(rem(it,2)~=0)%for sequential opertion bt EO and DE 
          
        for i = 1:nPop
               lambda=rand(1,Dim);                        
               r=rand(1,Dim);                             
               C_pool=[Ceq1; Ceq2; Ceq3; Ceq4; Ceq_ave];   % Equilibrium pool
               Ceq=C_pool(randi(size(C_pool,1)),:) ;       % random selection of one candidate from the pool
               F=a1*sign(r-0.5).*(exp(-lambda.*t)-1)  ;       
               r1=rand(); r2=rand();                          
               GCP=0.5*r1*ones(1,Dim)*(r2>=GP);                  
               G0=GCP.*(Ceq-lambda.*(pop(i).Position));
               G=G0.*F  ;                                       
            pop(i).Position = Ceq+(pop(i).Position-Ceq).*F+(G./lambda*V).*(1-F);
         end
    end   %end if 
   
  %% Show Iteration Information
   Convergence_curve(irun,it)=Ceq1_fit;
   end   %end it
 %% Show run Information
display(['Run no : ', num2str(irun)]);
display(['The best fitness found by EODE is : ', num2str(Ceq1_fit,10)]);
Bestfit_run(1,irun)=Ceq1_fit;
BestIndividual_run(irun,:)=Ceq1;

end   %end run
