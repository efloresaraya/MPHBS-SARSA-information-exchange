function [Best_Score, Best_Pos, Convergence_curve] = EODE_Wrapper(nPop, MaxIter, LB, UB, Dim, fobj)
    % Function wrapper for EODEsequential.m.
    % Adapted to accept an objective handle and dynamic bounds.
    
    % Internal EODE parameters.
    V=1; a1=2; a2=1; GP=0.5; FD=0.5; CR=0.9;
    
    % Initialization.
    DimSize = [1 Dim];
    Ceq1=zeros(1,Dim);   Ceq1_fit=inf; 
    Ceq2=zeros(1,Dim);   Ceq2_fit=inf; 
    Ceq3=zeros(1,Dim);   Ceq3_fit=inf; 
    Ceq4=zeros(1,Dim);   Ceq4_fit=inf;
    
    % Data structures.
    empty_individual.Position = [];
    empty_individual.Fitness = [];
    empty_individual.fit_old = [];
    empty_individual.Position_old = [];
    pop = repmat(empty_individual, nPop, 1);
    
    % Initialize population.
    for i = 1:nPop
        pop(i).Position = LB + rand(DimSize) .* (UB - LB);
        pop(i).Fitness = safe_eode_eval(fobj, pop(i).Position); 
    end
    
    Convergence_curve = zeros(1, MaxIter);
    
    % Main loop.
    for it = 1:MaxIter
        % Bound Constraint
        for i=1:nPop  
            Flag4ub=pop(i).Position>UB;
            Flag4lb=pop(i).Position<LB;
            pop(i).Position=(pop(i).Position).*(~(Flag4ub+Flag4lb))+UB.*Flag4ub+LB.*Flag4lb; 
        end
        
        % Evaluate fitness
        for i=1:nPop
            pop(i).Fitness = safe_eode_eval(fobj, pop(i).Position);
        end
        
        % Update Equilibrium Pool
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
        
        % Memory saving
        if it==1
          for i=1:nPop
             pop(i).fit_old=pop(i).Fitness;  pop(i).Position_old=pop(i).Position;
          end
        end
        for i=1:nPop
             if pop(i).fit_old < pop(i).Fitness
                 pop(i).Fitness=pop(i).fit_old; pop(i).Position=pop(i).Position_old;
             end
             pop(i).Position_old=pop(i).Position;  pop(i).fit_old=pop(i).Fitness;
        end
        
        % EODE equations.
        Ceq_ave=(Ceq1+Ceq2+Ceq3+Ceq4)/4;
        C_pool=[Ceq1; Ceq2; Ceq3; Ceq4; Ceq_ave];
        t=(1-it/MaxIter)^(a2*it/MaxIter);
        
        if (rem(it,2)==0) % DE Phase
           for i = 1:nPop
               rv1 = randi(nPop); while (rv1==i), rv1 = randi(nPop); end
               rv2 = randi(nPop); while ((rv2==rv1)||(rv2==i)), rv2 = randi(nPop); end
               rv3 = randi(nPop); while ((rv3==rv2)||(rv3==rv1)||(rv3==i)), rv3 = randi(nPop); end                  
               
               nrv=pop(rv3).Position+FD*(pop(rv1).Position-pop(rv2).Position);
               nrv = max(nrv, LB); nrv = min(nrv, UB); % Bounds check interno
               
               rndnum = rand(1,Dim) < CR;
               invmask = rndnum < 0.5;
               Trialv = invmask.*pop(i).Position + rndnum.*nrv;
               pop(i).Position=Trialv;
           end
        else % EO Phase
            for i = 1:nPop
               lambda=rand(1,Dim);                        
               r=rand(1,Dim);                             
               Ceq=C_pool(randi(size(C_pool,1)),:) ;       
               F=a1*sign(r-0.5).*(exp(-lambda.*t)-1);       
               r1=rand(); r2=rand();                          
               GCP=0.5*r1*ones(1,Dim)*(r2>=GP);                  
               G0=GCP.*(Ceq-lambda.*(pop(i).Position));
               G=G0.*F;                                       
               pop(i).Position = Ceq+(pop(i).Position-Ceq).*F+(G./lambda*V).*(1-F);
            end
        end 
        
        Convergence_curve(it) = Ceq1_fit;
    end
    
    Best_Score = Ceq1_fit;
    Best_Pos = Ceq1;
end

function y = safe_eode_eval(fobj, x)
    y = fobj(x(:).');
    if numel(y) > 1
        y = y(1);
    end
    if ~isfinite(y) || ~isreal(y)
        y = realmax;
    end
end
