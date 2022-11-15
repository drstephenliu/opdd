function checker = eval_tolerance(iter, objective, tolerance, range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% CHECK CONVERGENCE CRITERION, EXPRESSED AS THE AVERAGE RELATIVE CHANGE
%%%% OF OBJECTIVE FUNCTION OVER THE LAST N ITERATIONS (THEREFORE, AT LEAST 
%%%% N ITERATIONS ARE REQURED FOR TESTING CONVERGENCE
%%%% MAIN REFERENCE:
%%%%    [1] LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%% INPUT:
%%%%        ITER            --> CURRENT ITERATION COUNTER
%%%%        OBJECTIVE       --> OBJECTIVE FUNCTIONS OVER ITERATIONS
%%%%        TOLERANCE       --> TOLERANCE (BETWEEN 0 TO 100%)
%%%% OUTPUT:
%%%%        CHECKER         --> LOGICAL NUMBER INDICATING IF CONVERGED
%%%% CONTACT INFORMATION:
%%%%    STEPHEN Z. LIU 
%%%%    DEPARTMENT OF BIOMEDICAL ENGINEERING
%%%%    JOHNS HOPKINS UNIVERSITY SCHOOL OF MEDICINE
%%%%    E-MAIL: SZLIU29@GMAIL.COM
%%%%    PHONE: +1 657-275-2162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iter >= range 
    current_tol = mean(objective(iter-(range/2-1) : iter)); % CURRENT OBJECTIVE FUNCTION VALUE
    past_tol = mean(objective(iter-(range-1) : iter-(range/2))); % PREVIOUS OBJECTIVE FUNCTION VALUE
    relative_tol = abs(past_tol - current_tol) / current_tol; % RELATIVE CHANGE
    
    if relative_tol <= tolerance; checker = true;
    else; checker = false;
    end    
else
    checker = false;
end