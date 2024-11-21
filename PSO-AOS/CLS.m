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

% Input parameters
% Beta_k: The initial state of the chaotic map.
% FE: The current function evaluation number.
% FEmax: The predefined total number of function evaluations in an independent run.
% UB: Upper bound vector for solution candidates.
% LB: Lower bound vector for solution candidates.
% X: Selected solution candidate.
% Output
% X_CLS: The solution candidate derived from the CLS technique.

function [X_CLS, Beta_k]= CLS(Beta_k, FE, FEmax, UB, LB, X)
    Beta_k=2.3*(Beta_k.^2).*sin(pi*Beta_k);
    Lambda=(FEmax-FE+1)/FEmax;
    X_CLS=(1-Lambda)*X+Lambda*(LB+Beta_k*(UB-LB));
end