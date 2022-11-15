function varargout = get_fidelity(pie, y, M, B, M1, Z, DoGradient, DoHessian, DoObjective)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% COMPUTE GRADIENT, HESSIAN AND OBJECTIVE OF DATA FIDELITY TERM
%%%% MAIN REFERENCE:
%%%%        LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%% DATA FIDELITY IS EXPRESSED AS:
%%%%        L(PIE) = 1/2 * L2_NORM{Y - B*EXP(-M*PIE)}_Z
%%%% INPUT:
%%%%        PIE         --> VARIABLE
%%%%        Y           --> MEASUREMENT
%%%%        M           --> MATRIX OF LINEAR ATTENUATIONS
%%%%        B           --> MATRIX OF SPECTRAL RESPONSES
%%%%        M1          --> PRECOMPUTED PROJECTION OF ONES
%%%%        Z           --> INVERSE COVARIANCE (USE EYE IF NO NOISE MODEL)
%%%%        DOGRADIENT  --> IF COMPUTING GRADIENT (LOGICAL)
%%%%        DOHESSIAN   --> IF COMPUTING HESSIAN (LOGICAL)
%%%%        DOOBJECTIVE --> IF COMPUTING OBJECTIVE FUNCTION (LOGICAL)
%%%% OUTPUT:
%%%%        MAX OF THREE VARIABLES, IN THE ORDER OF 1) GRADIENT, 2) HESSIAN
%%%%        AND 3) OBJECTIVE FUNCTION.
%%%% CONTACT INFORMATION:
%%%%    STEPHEN Z. LIU 
%%%%    DEPARTMENT OF BIOMEDICAL ENGINEERING
%%%%    JOHNS HOPKINS UNIVERSITY SCHOOL OF MEDICINE
%%%%    E-MAIL: SZLIU29@GMAIL.COM
%%%%    PHONE: +1 657-275-2162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OO = exp(-M * pie);

if logical(DoGradient) % COMPUTE GRADIENT OF DATA FIDELITY
    garbage = y - (B * OO);
    garbage = Z .* garbage;
    garbage = B' * garbage;
    garbage = OO .* garbage;
    garbage = M' * garbage;
else; garbage = "none";
end

if logical(DoHessian) % COMPUTE (APPROXIMATED) HESSIAN OF DATA FIDELITY
    hessian = B * (OO .* M1);
    hessian = Z .* hessian;
    hessian = B' * hessian;
    hessian = OO .* hessian;
    hessian = M' * hessian;
else; hessian = "none";
end

if logical(DoObjective)  % COMPUTE OBJECTIVE FUNCTION OF DATA FIDELITY
    objective = y - (B * OO);
    objective = sqrt(Z) .* objective;
    objective = 0.5 * norm(objective, 'fro')^2;
else; objective = "none";
end

varargout{1} = garbage;
varargout{2} = hessian;
varargout{3} = objective;



