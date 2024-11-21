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

%  Honey Badger Algorithm source code 
%  paper:
%     Hashim, Fatma A., Essam H. Houssein, Kashif Hussain, Mai S. %     Mabrouk, Walid Al-Atabany.

function [Xprey, Food_Score,CNVG, Convergence_curve,AVG_Mem,FirstP_D1_Mem,DIV_Mem,Iter_Mem,Convergence_curve_Mem,AVG,Position_Mem] = HBAAOS(FunIndex,dim,down,up,T,Run_no,N,Beta_k)
for run=1:Run_no
beta       = 6;     % the ability of HB to get the food  Eq.(4)
C       = 2;     %constant in Eq. (3)
vec_flag=[1,-1];
%initialization
FE = 0;
FEmax = T * N;
average_objective = zeros(1, T);
FirstP_D1 = zeros(1 , T);
position_history = zeros(N , T , dim );
t = 0;
X=initialization(N,dim,up,down);
%Evaluation
for i =1:N
fitness(1, i) = BenFunctions((X(i,:)),FunIndex,dim);
FE = FE + 1;
end
[GYbest, gbest] = min(fitness);
Xprey = X(gbest,:);
while FE < FEmax
    t = t + 1;
    alpha=C*exp(-t/T);   %density factor in Eq. (3)
    I=Intensity(N,Xprey,X); %intensity in Eq. (2)
    for i=1:N
        position_history(i , t , : ) = X(i,:);
        r =rand();
        F=vec_flag(floor(2*rand()+1));
        for j=1:1:dim
            di=((Xprey(j)-X(i,j)));
            if r<.5
                r3=rand;                r4=rand;                r5=rand;
                
                Xnew(i,j)=Xprey(j) +F*beta*I(i)* Xprey(j)+F*r3*alpha*(di)*abs(cos(2*pi*r4)*(1-cos(2*pi*r5)));
            else
                r7=rand;
                Xnew(i,j)=Xprey(j)+F*r7*alpha*di;
            end
        end
        FU=Xnew(i,:)>up;FL=Xnew(i,:)<down;Xnew(i,:)=(Xnew(i,:).*(~(FU+FL)))+up.*FU+down.*FL;
        
        tempFitness = BenFunctions((Xnew(i,:)),FunIndex,dim);
        FE = FE + 1;
        if tempFitness<fitness(i)
            fitness(i)=tempFitness;
            X(i,:)= Xnew(i,:);
        end
           if FE >= FEmax
                 break
           end
    end
    FU=X>up;FL=X<down;X=(X.*(~(FU+FL)))+up.*FU+down.*FL;
    [Ybest,index] = min(fitness);
    CNVG(t)=min(Ybest);
    if Ybest<GYbest
        GYbest=Ybest;
        Xprey = X(index,:);
    end
    FirstP_D1(t) = X(1,1);
                      %-----------------------------------------------------------------------------------------------------------
[X, Xprey, GYbest, FE, Beta_k]= AOS(0.01, up.*ones(1,dim), down.*ones(1,dim), X,...
    fitness, FunIndex, dim, N, Xprey, GYbest, 1, dim, 1.5, Beta_k, FE, FEmax);
     for j=1:dim
         Div(j)=sum(median(X(:,j))-X(:,j))/N;
     end
     DIV(t)=mean(abs(Div));
       
       Convergence_curve(t)=GYbest; 
       for i=1:N
           average_objective(t) =  average_objective(t)  + fitness(i);
       end
       average_objective(t) = average_objective(t) / N;
           if FE >= FEmax
                 break
           end
    %-----------------------------------------------------------------------------------------------------------
end
Food_Score = GYbest;
res(run)=GYbest;
AVG_Mem(run)={average_objective};
FirstP_D1_Mem(run)={FirstP_D1};
Convergence_curve_Mem(run)={Convergence_curve};
Iter_Mem(run)={t-1};
DIV_Mem(run)={DIV};
Position_Mem(run)={position_history};
end
AVG=mean(res);
 fprintf('AVG F%s= %.1e\n', num2str(FunIndex), AVG);
end

function I=Intensity(N,Xprey,X)
for i=1:N-1
    di(i) =( norm((X(i,:)-Xprey+eps))).^2;
    S(i)=( norm((X(i,:)-X(i+1,:)+eps))).^2;
end
di(N)=( norm((X(N,:)-Xprey+eps))).^2;
S(N)=( norm((X(N,:)-X(1,:)+eps))).^2;
for i=1:N
    r2=rand;
    I(i)=r2*S(i)/(4*pi*di(i));
end
end
function [X]=initialization(N,dim,up,down)
if size(up,2)==1
    X=rand(N,dim).*(up-down)+down;
end
if size(up,2)>1
    for i=1:dim
        high=up(i);low=down(i);
        X(:,i)=rand(N,1).*(high-low)+low;
    end
end
end
