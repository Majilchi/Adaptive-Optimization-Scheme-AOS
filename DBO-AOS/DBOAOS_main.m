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
pop=30; % Number of search agents
M=500; % Maximum numbef of iterations
Run_no=30;
Beta_k=0.7;
func_num=50;
fun_to_B_run= [43, 43];
Conv=zeros(size(fun_to_B_run,2),Run_no,M);
Iter_M=zeros(size(fun_to_B_run,2),Run_no);
Min_Iter_Mem=zeros(1,size(fun_to_B_run,2));
% Load details of the selected benchmark functions
i=0;
for FunIndex = fun_to_B_run
i=i+1;
    [down,up,dim]=FunRange(FunIndex);
% for run = 1:30
[fMin,bestX,DBO_curve,AVG_Mem,FirstP_D1_Mem,DIV_Mem,Iter_Mem,Convergence_curve_Mem,AVG,Position_Mem]=DBOAOS(pop,M,down,up,dim,FunIndex,Run_no,Beta_k);
% res(run,1) = fMin;
% end
MemDBOAOS.AVG_Mem(i)={AVG_Mem};
MemDBOAOS.FirstP_D1_Mem(i)={FirstP_D1_Mem};
MemDBOAOS.DIV_Mem(i)={DIV_Mem};
MemDBOAOS.Iter_Mem(i)={Iter_Mem};
MemDBOAOS.Convergence_curve_Mem(i)={Convergence_curve_Mem};
MemDBOAOS.AVG(i)=AVG;
MemDBOAOS.Position_Mem(i)={Position_Mem};
end
for i=1:size(fun_to_B_run,2)
for j=1:Run_no
    for k=1:MemDBOAOS.Iter_Mem{1,i}{j}
        Conv(i,j,k)=MemDBOAOS.Convergence_curve_Mem{1,i}{1,j}(k);
    end
    Iter_M(i,j)=MemDBOAOS.Iter_Mem{1,i}{j};
end
figure
Min_Iter_Mem(1,i)=min(Iter_M(i,:));
plot(1:25:Min_Iter_Mem(1,i),squeeze(mean(Conv(i,:,1:25:Min_Iter_Mem(1,i)),2)),'-o','LineWidth',1,'MarkerFaceColor','auto','MarkerSize',5)
hold on
end