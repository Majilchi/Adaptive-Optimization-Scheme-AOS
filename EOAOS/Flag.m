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
% UB: Upper bound vector for solution candidates.
% LB: Lower bound vector for solution candidates.
% X: Selected solution candidate.
% dim: The number of dimensions of the selected solution candidate.
% Output
% X: The corrected solution candidate.

function X= Flag(UB, LB, X, dim)
for i=1:dim
    if (X(1,i)>UB(i))
        X(1,i)=UB(i);
    end
    if (X(1,i)<LB(i))
        X(1,i)=LB(i);
    end
end
end