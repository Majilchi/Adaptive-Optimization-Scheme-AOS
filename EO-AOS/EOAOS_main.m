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

clear all
clc
Run_no=30;         % Number of independent runs 
Particles_no=30;   % Number of particles
Max_iter=500; % Maximum number of iterations
Beta_k=0.7;
func_num=50;
fun_to_B_run= [43, 44];
Conv=zeros(size(fun_to_B_run,2),Run_no,Max_iter);
Iter_M=zeros(size(fun_to_B_run,2),Run_no);
Min_Iter_Mem=zeros(1,size(fun_to_B_run,2));
i=0;
for FunIndex= fun_to_B_run
i=i+1;      
[down,up,dim]=FunRange(FunIndex);
    
[Convergence_curve,Ave,Sd,AVG_Mem,FirstP_D1_Mem,DIV_Mem,Iter_Mem,Convergence_curve_Mem,Position_Mem]=EOAOS(Particles_no,Max_iter,down,up,dim,Run_no,FunIndex,Beta_k);
AVG(1,1)=Ave;
MemEOAOS.AVG_Mem(i)={AVG_Mem};
MemEOAOS.FirstP_D1_Mem(i)={FirstP_D1_Mem};
MemEOAOS.DIV_Mem(i)={DIV_Mem};
MemEOAOS.Iter_Mem(i)={Iter_Mem};
MemEOAOS.Convergence_curve_Mem(i)={Convergence_curve_Mem};
MemEOAOS.AVG(i)=AVG;
MemEOAOS.Position_Mem(i)={Position_Mem};
fprintf('AVG F%s= %.1e\n', num2str(FunIndex), AVG);
end
for i=1:size(fun_to_B_run,2)
for j=1:Run_no
    for k=1:MemEOAOS.Iter_Mem{1,i}{j}
        Conv(i,j,k)=MemEOAOS.Convergence_curve_Mem{1,i}{1,j}(k);
    end
    Iter_M(i,j)=MemEOAOS.Iter_Mem{1,i}{j};
end
figure
Min_Iter_Mem(1,i)=min(Iter_M(i,:));
plot(1:25:Min_Iter_Mem(1,i),squeeze(mean(Conv(i,:,1:25:Min_Iter_Mem(1,i)),2)),'-o','LineWidth',1,'MarkerFaceColor','auto','MarkerSize',5)
hold on
end