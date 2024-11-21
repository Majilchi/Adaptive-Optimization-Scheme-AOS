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

%  Equilibrium Optimizer source code (Developed in MATLAB R2015a)
%  programming: Afshin Faramarzi & Seyedali Mirjalili

function [Convergence_curve,Ave,Sd,AVG_Mem,FirstP_D1_Mem,DIV_Mem,Iter_Mem,Convergence_curve_Mem,Position_Mem]=EOAOS(Particles_no,Max_iter,down,up,dim,Run_no,FunIndex,Beta_k)
for irun=1:Run_no
Ceq1=zeros(1,dim);   Ceq1_fit=inf; 
Ceq2=zeros(1,dim);   Ceq2_fit=inf; 
Ceq3=zeros(1,dim);   Ceq3_fit=inf; 
Ceq4=zeros(1,dim);   Ceq4_fit=inf;
C=initialization(Particles_no,dim,up,down);
Iter=1; V=1;
a1=2;
a2=1;
GP=0.5;
FE=0;
FEmax=Max_iter*Particles_no;
average_objective = zeros(1, Max_iter);
FirstP_D1 = zeros(1 , Max_iter);
position_history = zeros(Particles_no , Max_iter , dim );
while FE<FEmax
      for i=1:size(C,1)  
        
        Flag4ub=C(i,:)>up;
        Flag4lb=C(i,:)<down;
        C(i,:)=(C(i,:).*(~(Flag4ub+Flag4lb)))+up.*Flag4ub+down.*Flag4lb;         
          position_history(i , Iter , : ) = C(i,:);
        fitness(i)=BenFunctions((C(i,:)),FunIndex,dim);
      
        if fitness(i)<Ceq1_fit
              Ceq1_fit=fitness(i);  Ceq1=C(i,:);
        elseif fitness(i)>Ceq1_fit && fitness(i)<Ceq2_fit  
              Ceq2_fit=fitness(i);  Ceq2=C(i,:);
        elseif fitness(i)>Ceq1_fit && fitness(i)>Ceq2_fit && fitness(i)<Ceq3_fit
              Ceq3_fit=fitness(i);  Ceq3=C(i,:);
        elseif fitness(i)>Ceq1_fit && fitness(i)>Ceq2_fit && fitness(i)>Ceq3_fit && fitness(i)<Ceq4_fit
              Ceq4_fit=fitness(i);  Ceq4=C(i,:);
                         
        end
        Ceqfit_run(irun)=Ceq1_fit;
        FE = FE + 1;
             if FE >= FEmax
                 break
             end
      end
       FirstP_D1(Iter) = C(1,1);
 %%%----AOS----%%%       
if (Iter>1)
    [C, Ceqfit_run, Ceq1, Ceq2, Ceq3, Ceq4, Ceq1_fit, Ceq2_fit, Ceq3_fit, Ceq4_fit, FE, Beta_k]=...
        AOS(0.01, up.*ones(1,dim), down.*ones(1,dim), C, fitness, FunIndex, dim, Particles_no, irun,...
        Ceqfit_run, Ceq1, Ceq2, Ceq3, Ceq4, Ceq1_fit, Ceq2_fit, Ceq3_fit, Ceq4_fit, 1, dim, 1.5, Beta_k, FE, FEmax);
end
      
%---------------- Memory saving-------------------   
      if Iter==1
        fit_old=fitness;  C_old=C;
      end
    
     for i=1:Particles_no
         if fit_old(i)<fitness(i)
             fitness(i)=fit_old(i); C(i,:)=C_old(i,:);
         end
     end
    C_old=C;  fit_old=fitness;
%-------------------------------------------------
       
Ceq_ave=(Ceq1+Ceq2+Ceq3+Ceq4)/4;                              % averaged candidate 
C_pool=[Ceq1; Ceq2; Ceq3; Ceq4; Ceq_ave];                     % Equilibrium pool
 
 t=(1-Iter/Max_iter)^(a2*Iter/Max_iter);                      % Eq (9)
 
    for i=1:Particles_no
           lambda=rand(1,dim);                                % lambda in Eq(11)
           r=rand(1,dim);                                     % r in Eq(11)  
           Ceq=C_pool(randi(size(C_pool,1)),:);               % random selection of one candidate from the pool
           F=a1*sign(r-0.5).*(exp(-lambda.*t)-1);             % Eq(11)
           r1=rand(); r2=rand();                              % r1 and r2 in Eq(15)
           GCP=0.5*r1*ones(1,dim)*(r2>=GP);                   % Eq(15)
           G0=GCP.*(Ceq-lambda.*C(i,:));                      % Eq(14)
           G=G0.*F;                                           % Eq(13)
           C(i,:)=Ceq+(C(i,:)-Ceq).*F+(G./lambda*V).*(1-F);   % Eq(16)                                                             
    end
     for j=1:dim
         Div(j)=sum(median(C(:,j))-C(:,j))/Particles_no;
     end
     DIV(Iter)=mean(abs(Div));
       
       Convergence_curve(Iter)=Ceq1_fit; 
       Ceqfit_run(irun)=Ceq1_fit;
       for i=1:Particles_no
           average_objective(Iter) =  average_objective(Iter)  + fitness(i);
       end
       average_objective(Iter) = average_objective(Iter) / Particles_no;
       Iter=Iter+1;  
end
AVG_Mem(irun)={average_objective};
FirstP_D1_Mem(irun)={FirstP_D1};
Convergence_curve_Mem(irun)={Convergence_curve};
Iter_Mem(irun)={Iter-1};
DIV_Mem(irun)={DIV};
Position_Mem(irun)={position_history};
end
Ave=mean(Ceqfit_run);
Sd=std(Ceqfit_run);  