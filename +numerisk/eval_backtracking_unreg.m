function alphas = eval_backtracking_unreg(pie, y, M, B, M1, Z, garbage, ...
    hessian, objective, Unit0, WolfeC1, WolfeC2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% BACKTRACKING LINE SEARCH USING NEWTON-RAPHSON (NR) SEARCH DIRECTION 
%%%% AND WOLFE CONDITIONS -- WITHOUT REGULARIZATION PART
%%%% MAIN REFERENCE:
%%%%        LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%% INPUT:
%%%%        PIE         --> VARIABLE
%%%%        Y           --> MEASUREMENT
%%%%        M           --> MATRIX OF LINEAR ATTENUATIONS
%%%%        B           --> MATRIX OF SPECTRAL RESPONSES
%%%%        M1          --> PRECOMPUTED PROJECTION OF ONES
%%%%        Z           --> INVERSE COVARIANCE (USE EYE IF NO NOISE MODEL)
%%%%        GARBAGE     --> GRADIENT AT CURRENT ITERATION
%%%%        HESSIAN     --> HESSIAN AT CURRENT ITERATION
%%%%        OBJECTIVE   --> OBJECTIVE FUNCTION AT CURRENT ITERATION
%%%%        UNIT0       --> INITIALIZED SEARCHING STEP
%%%%        WOLFEC1     --> CONSTANT FOR 1ST WOLFE CONDITION (CLOSE TO 0)
%%%%        WOLFEC2     --> CONSTANT FOR 2ND WOLFE CONDITION (CLOSE TO 1)
%%%% OUTPUT:
%%%%        ALPHAS      --> OPTIMAL STEPLENGTH
%%%% CONTACT INFORMATION:
%%%%    STEPHEN Z. LIU 
%%%%    DEPARTMENT OF BIOMEDICAL ENGINEERING
%%%%    JOHNS HOPKINS UNIVERSITY SCHOOL OF MEDICINE
%%%%    E-MAIL: SZLIU29@GMAIL.COM
%%%%    PHONE: +1 657-275-2162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE LINE SEARCH PARAMETERS
alpha0 = 1;
Unit0 = min(max(Unit0, 0), 1);
deltas = -1 * (1 ./ hessian) .* garbage; % NR SEARCH DIRECTION

% START BACKTRACKING
while alpha0 < inf
    pie_imagine = pie + alpha0 * deltas;
    
    % COMPUTE GRADIENT & OBJECTIVE OF TESTED UPDATE
    [garbage_imagine, ~, objective_imagine] = matematik.get_fidelity(pie_imagine, ...
                                                                     y, ...
                                                                     M, ...
                                                                     B, ...
                                                                     M1, ...
                                                                     Z, ...
                                                                     true, ...
                                                                     false, ...
                                                                     true);
    
    % EVALUATE WOLFE INEQUALITY
    orth = sum(deltas .* garbage, 'all');
    orth_imagine = sum(deltas .* garbage_imagine, 'all');
    wolfe1 = objective_imagine - objective;
    wolfe1 = wolfe1 - WolfeC1 * alpha0 * orth; % 1ST CONDITION (ARMIJO SUFFICIENT DECREASE)
    wolfe2 = abs(orth_imagine) - WolfeC2 * abs(orth); % 2ND CONDITION (WOLFE STRONG CURVATURE)
    
    % PASS THROUGH FILTER
    if (wolfe1 <= 0) && (wolfe2 <= 0) % STOP IFF BOTH CONDITIONS ARE SATISFIED
        alphas = alpha0; 
        break;
    else
        if (alpha0 - Unit0) <= 0; Unit0 = Unit0 .* 0.1; end % IF UNIT TOO LARGE, DECREASE BY 10%
        if Unit0 < 1e-4; alphas = 1e-5; break; end % IF UNIT TOO SMALL, TRUNCATE TO 0.00001
        alpha0 = alpha0 - Unit0; % OTHERWISE KEEP DECREASING STEPLENGTH
    end
end
