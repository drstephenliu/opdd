function varargout = optimize_linesearch_unreg(pie, y, M, B, Z, pPDD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% IMPLEMENT FAST NEWTON-RAPHSON OPTIMIZATION WITH LINE SEARCH & MOMENTUM
%%%% UPDATE
%%%% (NOTE THAT THIS IS UN-REGULARIZED VERSION. THE REGULARIZED VERSION IS
%%%% A DIFFERENT FILE)
%%%% MAIN REFERENCE:
%%%%    [1] LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%% INPUT:
%%%%        PIE             --> INITIALIZED VARIABLE 
%%%%        Y               --> MEASUREMENT
%%%%        M               --> MATRIX OF LINEAR ATTENUATIONS
%%%%        B               --> MATRIX OF SPECTRAL RESPONSES
%%%%        Z               --> MATRIX OF COVARIANCE
%%%%        PPDD            --> CLASS OF OPDD PARAMETERS
%%%% CONTACT INFORMATION:
%%%%    STEPHEN Z. LIU 
%%%%    DEPARTMENT OF BIOMEDICAL ENGINEERING
%%%%    JOHNS HOPKINS UNIVERSITY SCHOOL OF MEDICINE
%%%%    E-MAIL: SZLIU29@GMAIL.COM
%%%%    PHONE: +1 657-275-2162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHECK SOME USEFUL INPUTS
warning('on', 'all');
pPDD = numerisk.init_entries(pPDD);

% INITIALIZATION
fprintf('OPDD ->> INITIALIZING...\n');
objective_all = zeros(1, pPDD.ConvergeMaxIter); % ALLOCATE SPACE FOR OBJECTIVE
M1 = M * ones(size(pie)); % PRECOMPUTE PROJECTION OF ALL ONES
[pie0, time, time_sum, wind] = numerisk.init_momentum(pie, pPDD.UseGPU); % INITIALIZE MOMENTUM
if logical(pPDD.UseGPU); objective_all = gpuArray(objective_all); M1 = gpuArray(M1); end

% START OPTIMIZATION
fprintf('OPDD ->> START OPTIMIZING...\n');
for iter = 1 : pPDD.ConvergeMaxIter
    if rem(iter, 1000) == 0; disp(['ITER = ' num2str(iter)]); end
    
    % COMPUTE DATA FIDELITY GRADIENT AND HESSIAN
    [garbage, hessian, objective] = matematik.get_fidelity(pie, ...
                                                           y, ...
                                                           M, ...
                                                           B, ...
                                                           M1, ...
                                                           Z, ...
                                                           true, ...
                                                           true, ...
                                                           true);
    
    objective_all(:, iter) = objective; % SAVE OBJECTIVE FOR CONVERGENCE CHECKING
    
    % LINE SEARCH (ONE STEPLENGTH FOR ALL PIXELS)
    alphas = numerisk.eval_backtracking_unreg(pie, ...
                                              y, ...
                                              M, ...
                                              B, ...
                                              M1, ...
                                              Z, ...
                                              garbage, ...
                                              hessian, ...
                                              objective, ...
                                              pPDD.BacktrackInitStep, ...
                                              pPDD.BacktrackWolfeC1, ...
                                              pPDD.BacktrackWolfeC2);
    
    % MOMENTUM/REGULAR UPDATE
    deltas = alphas * (1 ./ hessian) .* garbage;
    [pie, pie0, time, time_sum, wind] = numerisk.eval_momentum(iter, ...
                                                               deltas, ...
                                                               pie, ...
                                                               pie0, ...
                                                               time, ...
                                                               time_sum, ...
                                                               wind, ...
                                                               pPDD.MomentMemory, ...
                                                               pPDD.MomentTruncate(1), ...
                                                               pPDD.MomentTruncate(2));
    
    % CHECKING TOLERANCE
    if numerisk.eval_tolerance(iter, objective_all, pPDD.ConvergeTolerance, pPDD.ConvergeRange); break; 
    end
    
    % PLOT CONVERGENCE IF YOU LIKE
    if logical(pPDD.DrawConvergence)
        if rem(iter, pPDD.DrawUpdateIter) == 1
            loglog(objective_all(1 : iter)', 'LineWidth', 1.5);
            legend('OBJECTIVE FUNCTION');
            title(['CONVERGENCE --> ITERATION = ' num2str(iter)]);
            xlabel('ITERATION'); ylabel('OBJECTIVE FUNCTION'); 
            set(gca, 'FontSize', 14, 'FontName', 'Courier');
            drawnow;
        end
    end
end
drawnow;
varargout{1} = pie;
varargout{2} = objective_all;
varargout{3} = pPDD;