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

% Dung Beetle Optimizer: (DBO) (demo)
% Programmed by Jian-kai Xue

function [fMin , bestX, Convergence_curve,AVG_Mem,FirstP_D1_Mem,DIV_Mem,Iter_Mem,Convergence_curve_Mem,AVG,Position_Mem ] = DBOAOS(pop,M,down,up,dim,FunIndex,Run_no,Beta_k)
for run=1:Run_no        
   P_percent = 0.2;    % The population size of producers accounts for "P_percent" percent of the total population size       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pNum = round( pop *  P_percent );    % The population size of the producers   
lb= down.*ones( 1,dim );    % Lower limit/bounds/     a vector
ub= up.*ones( 1,dim );    % Upper limit/bounds/     a vector
%----------------------------------------------------------------
FE = 0;
FEmax = pop * M;
t = 0;
average_objective = zeros(1, M);
FirstP_D1 = zeros(1 , M);
position_history = zeros(pop , M , dim );
%----------------------------------------------------------------
%Initialization
for i = 1 : pop
    x( i, : ) = lb + (ub - lb) .* rand( 1, dim );  
    fit( i ) = BenFunctions( x( i, : ),FunIndex,dim ) ;     
    FE = FE + 1;
end
pFit = fit;                       
pX = x; 
 XX=pX;    
[ fMin, bestI ] = min( fit );      % fMin denotes the global optimum fitness value
bestX = x( bestI, : );             % bestX denotes the global optimum position corresponding to fMin
 % Start updating the solutions.
while FE<FEmax 
         t = t + 1;
        [fmax,B]=max(fit);
        worse= x(B,:);   
       r2=rand(1);
       for i=1:pop
           position_history(i , t , : ) = x(i,:);
       end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : pNum    
        if(r2<0.9)
            r1=rand(1);
          a=rand(1,1);
          if (a>0.1)
           a=1;
          else
           a=-1;
          end
    x( i , : ) =  pX(  i , :)+0.3*abs(pX(i , : )-worse)+a*0.1*(XX( i , :)); % Equation (1)
       else
            
           aaa= randperm(180,1);
           if ( aaa==0 ||aaa==90 ||aaa==180 )
            x(  i , : ) = pX(  i , :);   
           end
         theta= aaa*pi/180;   
       
       x(  i , : ) = pX(  i , :)+tan(theta).*abs(pX(i , : )-XX( i , :));    % Equation (2)      
        end
      
        x(  i , : ) = Bounds( x(i , : ), lb, ub );    
        fit(  i  ) = BenFunctions( x( i, : ),FunIndex,dim );
        FE = FE +1;
    end 
 [ fMMin, bestII ] = min( fit );      % fMin denotes the current optimum fitness value
  bestXX = x( bestII, : );             % bestXX denotes the current optimum position 
 R=1-t/M;                           %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Xnew1 = bestXX.*(1-R); 
     Xnew2 =bestXX.*(1+R);                    %%% Equation (3)
   Xnew1= Bounds( Xnew1, lb, ub );
   Xnew2 = Bounds( Xnew2, lb, ub );
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Xnew11 = bestX.*(1-R); 
     Xnew22 =bestX.*(1+R);                     %%% Equation (5)
   Xnew11= Bounds( Xnew11, lb, ub );
    Xnew22 = Bounds( Xnew22, lb, ub );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i = ( pNum + 1 ) :12                  % Equation (4)
     x( i, : )=bestXX+((rand(1,dim)).*(pX( i , : )-Xnew1)+(rand(1,dim)).*(pX( i , : )-Xnew2));
   x(i, : ) = Bounds( x(i, : ), Xnew1, Xnew2 );
  fit(i ) = BenFunctions( x( i, : ),FunIndex,dim ) ;
  FE = FE +1;
   end
   
  for i = 13: 19                  % Equation (6)
   
        x( i, : )=pX( i , : )+((randn(1)).*(pX( i , : )-Xnew11)+((rand(1,dim)).*(pX( i , : )-Xnew22)));
       x(i, : ) = Bounds( x(i, : ),lb, ub);
       fit(i ) = BenFunctions( x( i, : ),FunIndex,dim ) ;
       FE = FE +1;
  end
  
  for j = 20 : pop                 % Equation (7)
       x( j,: )=bestX+randn(1,dim).*((abs(( pX(j,:  )-bestXX)))+(abs(( pX(j,:  )-bestX))))./2;
      x(j, : ) = Bounds( x(j, : ), lb, ub );
      fit(j ) = BenFunctions( x( j, : ),FunIndex,dim ) ;
      FE = FE +1;
  end
   % Update the individual's best fitness vlaue and the global best fitness value
     XX=pX;
    for i = 1 : pop 
        if ( fit( i ) < pFit( i ) )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end
        
        if( pFit( i ) < fMin )
           % fMin= pFit( i );
            fMin= pFit( i );
            bestX = pX( i, : );
          %  a(i)=fMin;
            
        end
    end
    FirstP_D1(t) = x(1,1);
                  %-----------------------------------------------------------------------------------------------------------
[x, bestX, XX, pX, fMin, pFit, FE, Beta_k]=AOS(0.01, ub,...
    lb, x, fit, FunIndex, dim, pop, bestX,...
    XX, pX, fMin, pFit, 1, dim, 1.5, Beta_k, FE, FEmax);


    %-----------------------------------------------------------------------------------------------------------
  
     
     for j=1:dim
         Div(j)=sum(median(x(:,j))-x(:,j))/pop;
     end
     DIV(t)=mean(abs(Div));
       
       Convergence_curve(t)=fMin; 
       for i=1:pop
           average_objective(t) =  average_objective(t)  + fit(i);
       end
       average_objective(t) = average_objective(t) / pop;
    
           if FE >= FEmax
                 break
           end
end
res(run)=fMin;
AVG_Mem(run)={average_objective};
FirstP_D1_Mem(run)={FirstP_D1};
Convergence_curve_Mem(run)={Convergence_curve};
Iter_Mem(run)={t-1};
DIV_Mem(run)={DIV};
Position_Mem(run)={position_history};
end

AVG=mean(res);
fprintf('AVG F%s= %.1e\n', num2str(FunIndex), AVG);
% Application of simple limits/bounds
function s = Bounds( s, Lb, Ub)
  % Apply the lower bound vector
  temp = s;
  I = temp < Lb;
  temp(I) = Lb(I);
  
  % Apply the upper bound vector 
  J = temp > Ub;
  temp(J) = Ub(J);
  % Update this new move 
  s = temp;
function S = Boundss( SS, LLb, UUb)
  % Apply the lower bound vector
  temp = SS;
  I = temp < LLb;
  temp(I) = LLb(I);
  
  % Apply the upper bound vector 
  J = temp > UUb;
  temp(J) = UUb(J);
  % Update this new move 
  S = temp;
%---------------------------------------------------------------------------------------------------------------------------
