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

clear 
close all
clc
noP = 30;
T = 500;
Run_no=30;
BestSolutions_PSO = zeros(1 , Run_no);
Beta_k=0.7;
func_num=50;
fun_to_B_run= [43];
Conv=zeros(size(fun_to_B_run,2),Run_no,T);
Iter_M=zeros(size(fun_to_B_run,2),Run_no);
Min_Iter_Mem=zeros(1,size(fun_to_B_run,2));
i=0;
for FunIndex= fun_to_B_run
    i=i+1;
    [down,up,dim]=FunRange(FunIndex);
    [ GBEST, Convergence_curve,AVG_Mem,FirstP_D1_Mem,DIV_Mem,Iter_Mem,Convergence_curve_Mem,AVG,Position_Mem] = PSOAOS( noP , T, dim, up, down, FunIndex, Run_no, Beta_k) ;
MemPSOAOS.AVG_Mem(i)={AVG_Mem};
MemPSOAOS.FirstP_D1_Mem(i)={FirstP_D1_Mem};
MemPSOAOS.DIV_Mem(i)={DIV_Mem};
MemPSOAOS.Iter_Mem(i)={Iter_Mem};
MemPSOAOS.Convergence_curve_Mem(i)={Convergence_curve_Mem};
MemPSOAOS.AVG(i)=AVG;
MemPSOAOS.Position_Mem(i)={Position_Mem};
end
for i=1:size(fun_to_B_run,2)
for j=1:Run_no
    for k=1:MemPSOAOS.Iter_Mem{1,i}{j}
        Conv(i,j,k)=MemPSOAOS.Convergence_curve_Mem{1,i}{1,j}(k);
    end
    Iter_M(i,j)=MemPSOAOS.Iter_Mem{1,i}{j};
end
figure
Min_Iter_Mem(1,i)=min(Iter_M(i,:));
plot(1:25:Min_Iter_Mem(1,i),squeeze(mean(Conv(i,:,1:25:Min_Iter_Mem(1,i)),2)),'-o','LineWidth',1,'MarkerFaceColor','auto','MarkerSize',5)
hold on
end