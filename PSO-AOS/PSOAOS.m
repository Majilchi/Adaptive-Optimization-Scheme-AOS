%--------------------------------------------------------------------------------------------------------------------------
% Adaptive optimization scheme (AOS)
% Majid Ilchi Ghazaan, Amirmohammad Salmani Oshnari, Amirhossein Salmani Oshnari,
% A novel adaptive optimization scheme for advancing metaheuristics and global optimization,
% Swarm and Evolutionary Computation, Volume 91, 2024, 101779, ISSN 2210-6502,
% DOI: https://doi.org/10.1016/j.swevo.2024.101779.
% (https://www.sciencedirect.com/science/article/pii/S2210650224003171)
% Programmed by: Amirmohammad Salmani Oshnari & Amirhossein Salmani Oshnari
% e-mail: ilchi@iust.ac.ir, amirmohammad78.salmani@gmail.com & salmaniamirhossein78@gmail.com
% Updated 21 Nov. 2024.
%--------------------------------------------------------------------------------------------------------------------------

% Written by Dr. Seyedali Mirjalili


function [ GBEST, Convergence_curve,AVG_Mem,FirstP_D1_Mem,DIV_Mem,Iter_Mem,Convergence_curve_Mem,AVG,Position_Mem] = PSOAOS( noP , T , dim, up, down, FunIndex, Run_no, Beta_k)
for run=1:Run_no
lb=down.*ones(1,dim);
ub=up.*ones(1,dim);
average_objective = zeros(1, T);
FirstP_D1 = zeros(1 , T);
position_history = zeros(noP , T , dim );
SwarmPX= zeros(noP, dim);
SwarmPO= zeros(noP, 1);
wMax = 0.9;
wMin = 0.2;
c1 = 2;
c2 = 2;
vMax = (ub - lb) .* 0.2;
vMin  = -vMax;
FE=0;
FEmax=T*noP;
t=0;
% The PSO algorithm

% Initialize the particles
for k = 1 : noP
    Swarm.Particles(k).X = (ub-lb) .* rand(1,dim) + lb;
    Swarm.Particles(k).V = zeros(1, dim);
    Swarm.Particles(k).PBEST.X = zeros(1,dim);
    Swarm.Particles(k).PBEST.O = inf;
    
    Swarm.GBEST.X = zeros(1,dim);
    Swarm.GBEST.O = inf;
end


% Main loop
while(FE<FEmax)
    t=t+1;
    % Calcualte the objective value
    for k = 1 : noP
        
        position_history(k , t , : ) = Swarm.Particles(k).X;
        
        Swarm.Particles(k).O = BenFunctions(Swarm.Particles(k).X, FunIndex, dim);
             FE = FE + 1;
         % Update the PBEST
        if Swarm.Particles(k).O < Swarm.Particles(k).PBEST.O
            Swarm.Particles(k).PBEST.X = Swarm.Particles(k).X;
            Swarm.Particles(k).PBEST.O = Swarm.Particles(k).O;
        end
        
        % Update the GBEST
        if Swarm.Particles(k).O < Swarm.GBEST.O
            Swarm.GBEST.X = Swarm.Particles(k).X;
            Swarm.GBEST.O = Swarm.Particles(k).O;
        end
        GBEST = Swarm.GBEST;
             if FE >= FEmax
                 break
             end
             SwarmPX(k, :)= Swarm.Particles(k).X;
             SwarmPO(k,1)= Swarm.Particles(k).O;
    end
SwarmBestX= Swarm.GBEST.X;
SwarmBestO= Swarm.GBEST.O;
if t>1
    [SwarmPX, SwarmBestX, SwarmBestO, FE, Beta_k]=AOS(0.01, ub, lb, SwarmPX, SwarmPO,...
        FunIndex, dim, noP, SwarmBestX, SwarmBestO, 1, dim, 1.5, Beta_k, FE, FEmax);
end
Swarm.GBEST.X= SwarmBestX;
Swarm.GBEST.O= SwarmBestO;
    %-----------------------------------------------------------------------------------------------------------
    for k = 1 : noP
        Swarm.Particles(k).X= SwarmPX(k);
        Swarm.Particles(k).O= SwarmPO(k);
        % Update the PBEST
        if Swarm.Particles(k).O < Swarm.Particles(k).PBEST.O
            Swarm.Particles(k).PBEST.X = Swarm.Particles(k).X;
            Swarm.Particles(k).PBEST.O = Swarm.Particles(k).O;
        end
        
        % Update the GBEST
        if Swarm.Particles(k).O < Swarm.GBEST.O
            Swarm.GBEST.X = Swarm.Particles(k).X;
            Swarm.GBEST.O = Swarm.Particles(k).O;
        end
    end
    
    % Update the X and V vectors
    w = wMax - t .* ((wMax - wMin) / T);
    
    for k = 1 : noP
        Swarm.Particles(k).V = w .* Swarm.Particles(k).V + c1 .* rand(1,dim) .* (Swarm.Particles(k).PBEST.X - Swarm.Particles(k).X) ...
            + c2 .* rand(1,dim) .* (Swarm.GBEST.X - Swarm.Particles(k).X);
        
        
        % Check velocities
        index1 = find(Swarm.Particles(k).V > vMax);
        index2 = find(Swarm.Particles(k).V < vMin);
        
        Swarm.Particles(k).V(index1) = vMax(index1);
        Swarm.Particles(k).V(index2) = vMin(index2);
        
        Swarm.Particles(k).X = Swarm.Particles(k).X + Swarm.Particles(k).V;
        
        % Check positions
        index1 = find(Swarm.Particles(k).X > ub);
        index2 = find(Swarm.Particles(k).X < lb);
        
        Swarm.Particles(k).X(index1) = ub(index1);
        Swarm.Particles(k).X(index2) = lb(index2);
        
    end
        FirstP_D1(t) = Swarm.Particles(1).X(1);
    for i=1:noP
        X_Plot(i,:)=Swarm.Particles(i).X;
    end
     for j=1:dim
         Div(j)=sum(median(X_Plot(:,j))-X_Plot(:,j))/noP;
     end
     DIV(t)=mean(abs(Div));
       
       Convergence_curve(t)=Swarm.GBEST.O; 
       for i=1:noP
           average_objective(t) =  average_objective(t)  + Swarm.Particles(i).O;
       end
       average_objective(t) = average_objective(t) / noP;

end

GBEST = Swarm.GBEST;
res(run)=Swarm.GBEST.O;
AVG_Mem(run)={average_objective};
FirstP_D1_Mem(run)={FirstP_D1};
Convergence_curve_Mem(run)={Convergence_curve};
Iter_Mem(run)={t-1};
DIV_Mem(run)={DIV};
Position_Mem(run)={position_history};
end

AVG=mean(res);
fprintf('AVG F%s= %.1e\n', num2str(FunIndex), AVG);