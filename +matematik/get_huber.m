function varargout = get_huber(pie, HuberDelta, HuberBeta, Dim, DoGradient, DoHessian, DoObjective)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% COMPUTE GRADIENT, HESSIAN AND OBJECTIVE OF REGULARIZATION TERM
%%%% MAIN REFERENCE:
%%%%        LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%% DATA FIDELITY IS EXPRESSED AS:
%%%%        R(PIE) = SUM(BETA * HUBER(PIE))
%%%% INPUT:
%%%%        PIE         --> VARIABLE
%%%%        HUBERDELTA  --> HUBER CUTOFF (CURRENTLY A SCALAR)
%%%%        HUBERBETA   --> HUBER STRENGTH (SAME LENGTH AS #MATERIALS)
%%%%        DIM         --> SIZE OF RADIOGRAPH IN 2D
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

pie = reshape(pie, [], Dim(1), Dim(2)); % RESHAPE PIE TO #MATERIAL X 2D SIZE

% JUST TO ALLOCATE SPACE FOR OUTPUTS
garbage = pie .* 0; 
hessian = pie .* 0;
objective = 0;

% CHECK IF NUMBER OF HYPERPARAMETERS ARE ENOUGH
if ~isscalar(HuberDelta); error('HUBERDELTA NOT A SCALAR...'); end
if isscalar(HuberBeta); HuberBeta = repmat(HuberBeta, size(pie, 1), 1);
else; HuberBeta = reshape(HuberBeta, [], 1);
end

% SOME RANDOM FUNCTIONS FOR FINITE DIFFERENCE
link1 = @(O) padarray(O(:, 2 : end, :) - O(:, 1 : end-1, :), [0 1 0], 0, 'pre');
link2 = @(O) padarray(O(:, 1 : end-1, :) - O(:, 2 : end, :), [0 1 0], 0, 'post');
link3 = @(O) padarray(O(:, :, 2 : end) - O(:, :, 1 : end-1), [0 0 1], 0, 'pre');
link4 = @(O) padarray(O(:, :, 1 : end-1) - O(:, :, 2 : end), [0 0 1], 0, 'post');
if size(pie, 3) == 1; linker = {link1, link2}; else; linker = {link1, link2, link3, link4}; end
if size(pie, 3) == 1; kvarter = 2; else; kvarter = 4; end

% MAIN PART FOR COMPUTING HUBER REGULARIZATION TERMS
for kay = 1 : kvarter
    link_kay = linker{kay}(pie);
    cut_kay = abs(link_kay) < HuberDelta;
    
    if logical(DoGradient)
        corn = link_kay / HuberDelta;
        corn = corn .* cut_kay + sign(corn) .* ~cut_kay;
        garbage = garbage + corn;
    end
    
    if logical(DoHessian)
        corn = link_kay .* 0 + (1 ./ HuberDelta);
        corn = corn .* cut_kay;
        hessian = hessian + corn;
    end 
    
    if logical(DoObjective)
        corn = (link_kay.^2 ./ (2 * HuberDelta)) .* cut_kay;
        corn = corn + (abs(link_kay) - HuberDelta / 2) .* ~cut_kay;
        objective = objective + HuberBeta' * sum(corn, [2 3]);
    end 
end

% WRAP UP ALL OUTPUTS
if ~logical(DoGradient); garbage = "none"; 
else; garbage = reshape(garbage, [], Dim(1)*Dim(2)) .* HuberBeta;
end
if ~logical(DoHessian); hessian = "none"; 
else; hessian = reshape(hessian, [], Dim(1)*Dim(2)) .* HuberBeta;
end
if ~logical(DoObjective); objective = "none"; end

varargout{1} = garbage;
varargout{2} = hessian;
varargout{3} = objective;



