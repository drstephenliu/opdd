%% OPPD WITH OPTIMAL LINE SEARCH, MOMENTUM AND HUBER REGULARIZATION
clearvars
clc
close all

reset(gpuDevice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% EXAMPLE OF THREE-ENERGY TWO-MATERIAL DECOMPOSITION FOR THE DATA
%%%% OBTAINED FROM TRIPLE-LAYER FLAT-PANEL DETECTOR
%%%% HERE, WE LISTED ALL INPUTS.
%%%% HOWEVER, DEFAULT INPUT (SEE THE OTHER EXAMPLE SCRIPT) IS RECOMMENDED
%%%% TO START WITH.
%%%% (NOTE THAT UN-REGULARIZED VERSION IS A DIFFERENT FILE)
%%%% MAIN REFERENCE:
%%%%    [1] LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%% CONTACT INFORMATION:
%%%%    STEPHEN Z. LIU 
%%%%    DEPARTMENT OF BIOMEDICAL ENGINEERING
%%%%    JOHNS HOPKINS UNIVERSITY SCHOOL OF MEDICINE
%%%%    E-MAIL: SZLIU29@GMAIL.COM
%%%%    PHONE: +1 657-275-2162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DATA FILES
pFNS.Data1              = 'data/proj_layer1_120kVp_1500uAs_grid_heel.mat'; % 1ST LAYER DATA
pFNS.Data2              = 'data/proj_layer2_120kVp_1500uAs_grid_heel.mat'; % 2ND LAYER DATA
pFNS.Data3              = 'data/proj_layer3_120kVp_1500uAs_grid_heel.mat'; % 3RD LAYER DATA
pFNS.OPDD1              = 'demo/opdd_reg_water_radiograph_customize.mat'; % FOR SAVING, 1ST MATERIAL RADIOGRAPH 
pFNS.OPDD2              = 'demo/opdd_reg_bone_radiograph_customize.mat'; % FOR SAVING, 2ND MATERIAL RADIOGRAPH 
                       
pGEO.Source             = 'physics/spectrum_120kVp_2000al+250cu.mat'; % SOURCE SPECTRUM (1-150 KEV ARRAY)
pGEO.Grid               = 'physics/spectrum_grid_vertical_13to1_pb+al.mat'; % ANTISCATTER GRID SPECTRUM (1-150 KEV ARRAY)
pGEO.Detector1          = 'physics/detector_layer1_210csi.mat'; % 1ST LAYER DETECTOR RESPONSE (1-150 KEV ARRAY)
pGEO.Detector2          = 'physics/detector_layer2_210csi+560csi.mat'; % 2ND LAYER DETECTOR RESPONSE (1-150 KEV ARRAY)
pGEO.Detector3          = 'physics/detector_layer3_210csi+560csi+550csi.mat'; % 3RD LAYER DETECTOR RESPONSE (1-150 KEV ARRAY)
pGEO.Attenuation1       = 'physics/linear_attenuation_water.mat'; % LINEAR ATTENUATION OF 1ST MATERIAL (1-150 KEV ARRAY)
pGEO.Attenuation2       = 'physics/linear_attenuation_bone.mat'; % LINEAR ATTENUATION OF 2ND MATERIAL (1-150 KEV ARRAY)

% OPDD HYPERPARAMETERS (HERE, ALL INPUTS ARE LISTED)
pPDD.UseGPU             = true; % IF USING GPU
pPDD.UseMomentum        = true; % IF USING MOMENTUM UPDATE
pPDD.UseHuber           = true; % IF USING REGULARIZATION
pPDD.UseHuberAdjust     = true; % IF USING DYNAMIC REGULARIZATION STRENGTH ADJUSTMENT
pPDD.ConvergeTolerance  = 0.0001; % CONVERGENCE TOLERANCE (FRACTION OF CURRENT OBJECTIVE VALUE)
pPDD.ConvergeMaxIter    = 2000; % MAXIMUM ITERATIONS
pPDD.ConvergeRange      = 20; % NUMBER OF MOST RECENT ITERATIONS FOR COMPUTING CONVERGENCE CRITERION
pPDD.BacktrackInitStep  = 0.1; % LINE SEARCH INITIAL SEARCHING STEP (VALUE BETWEEN 0 TO 1)
pPDD.BacktrackWolfeC1   = 0.0001; % 1ST WOLFE (ARMIJO) CONDITION CONSTANT, USUALLY CLOSE TO 0
pPDD.BacktrackWolfeC2   = 0.9999; % 2ND WOLFE (STRONG) CONDITION CONSTANT, USUALLY CLOSE TO 1
pPDD.MomentMemory       = 100; % LENGTH THAT MOMENTUM MEMORIZES, RESET AFTER THIS
pPDD.MomentTruncate     = [-10 inf]; % LOWER AND UPPER TRUNCATION FOR VARIABLE (E.G., HERE, I ADDED A LOWER TRUNCATION OF -10)
pPDD.HuberImageSize     = [450 288]; % RADIOGRAPH SIZE IN PIXELS (NECESSARY INPUT IF USEHUBER)
pPDD.HuberBeta          = [1e-7 1e-6]; % INITIAL REGULARIZATION STRENGTH PER MATERIAL (NECESSARY INPUT IF USEHUBER)
pPDD.HuberDelta         = 0.5; % HUBER CUTOFF (NECESSARY INPUT IF USEHUBER)
pPDD.HuberAdjustRatio   = 0.4; % RATIO OF REGULARIZATION OVER TOTAL OBJECTIVE, FOR ADJUSTING HUBER STRENGTH
pPDD.HuberAdjustFreq    = 200; % NUMBER OF ITERATIONS TO ADJUST HUBER STRENGTH
pPDD.HuberAdjustRange   = 10; % NUMBER OF MOST RECENT ITERATIONS FOR COMPUTING HUBER ADJUST RATIO
pPDD.DrawConvergence    = true; % IF PLOTTING CONVERGENCE
pPDD.DrawUpdateIter     = 100; % NUMBER OF ITERATIONS TO UPDATE CONVERGENCE PLOT

% NOTE THAT HERE, SINCE WE USE HUBER ADJUST, THE INITIAL INPUT FOR BETA CAN
% BE ANY VALUES (PREFERABLY CLOSE TO ZERO). THIS IS BECAUSE THE VALUES WILL
% BE AUTOMATICALLY COMPUTED AND ADJUSTED AFTER CERTAIN ITERATIONS.

% =========================================================================
% ===================== NO INPUT STARTING FROM HERE =======================
% =========================================================================

% BUILD DATA VECTOR (Y VECTOR)
Y1 = struct2array(load(pFNS.Data1));
Y2 = struct2array(load(pFNS.Data2));
Y3 = struct2array(load(pFNS.Data3));
Y = gpuArray([reshape(Y1, [], 1)'; reshape(Y2, [], 1)'; reshape(Y3, [], 1)']);


% BUILDING SPECTRUM MATRIX (B MATRIX)
spec = struct2array(load(pGEO.Source));
grd = struct2array(load(pGEO.Grid));
det1 = struct2array(load(pGEO.Detector1));
det2 = struct2array(load(pGEO.Detector2));
det3 = struct2array(load(pGEO.Detector3));
s1 = det1 .* spec .* grd;
s2 = det2 .* spec .* grd;
s3 = det3 .* spec .* grd;
s1 = s1 ./ sum(s1);
s2 = s2 ./ sum(s2);
s3 = s3 ./ sum(s3);
s1 = s1(5:120); % REMOVE 1-4 KEV POINTS, WHICH ARE USELESS AND INACCURATE
s2 = s2(5:120); % REMOVE 1-4 KEV POINTS, WHICH ARE USELESS AND INACCURATE
s3 = s3(5:120); % REMOVE 1-4 KEV POINTS, WHICH ARE USELESS AND INACCURATE
B = gpuArray([s1'; s2'; s3']);


% BUILDING MATERIAL MATRIX (M MATRIX)
atten1 = struct2array(load(pGEO.Attenuation1));
atten2 = struct2array(load(pGEO.Attenuation2));
atten1 = atten1(5:120); % REMOVE 1-4 KEV POINTS, WHICH ARE USELESS AND INACCURATE
atten2 = atten2(5:120); % REMOVE 1-4 KEV POINTS, WHICH ARE USELESS AND INACCURATE
M = gpuArray([atten1 atten2]);


% BUILDING COVARIANCE WEIGHTING MATRIX (Z MATRIX)
Z = 1 ./ Y;


% INITIALIZING VARIABLE (PIE VECTOR)
pie = gpuArray(ones(2, size(Y, 2))); % HERE, FOR EXAMPLE, I INITIALIZE WITH ALL 1


% RUNNING OPTIMIZATION
[pie, objective_all, pPDD] = opdd.optimize_linesearch_reg(pie, Y, M, B, Z, pPDD);


% WRAP UP RESULTS
pie = gather(pie);
pie1 = reshape(pie(1, :), size(Y1, 1), size(Y1, 2)); % JUST TO RESHAPE FROM VECTOR BACK TO RADIOGRAPH
pie2 = reshape(pie(2, :), size(Y1, 1), size(Y1, 2)); % JUST TO RESHAPE FROM VECTOR BACK TO RADIOGRAPH
save(pFNS.OPDD1, 'pie1');
save(pFNS.OPDD2, 'pie2');
disp('DONE!');


%% PLOTTING RESULT
f1 = figure; imagesc(fliplr(rot90(pie1))); 
axis image; colormap(gray); caxis([0 250]); title('Decomposed water pathlength');
set(gca, 'FontSize', 16, 'FontName', 'Courier', 'LineWidth', 1.5, 'XTickLabel', {}, 'YTickLabel', {}, 'XColor', 'k', 'YColor', 'k');
c = colorbar(); c.LineWidth = 1.5; c.FontSize = 20; c.Location = 'southoutside'; c.Color = 'k';
saveas(f1, [pFNS.OPDD1(1 : end-4), '.png']);

f2 = figure; imagesc(fliplr(rot90(pie2))); 
axis image; colormap(gray); caxis([-10 10]); title('Decomposed bone pathlength');
set(gca, 'FontSize', 16, 'FontName', 'Courier', 'LineWidth', 1.5, 'XTickLabel', {}, 'YTickLabel', {}, 'XColor', 'k', 'YColor', 'k');
c = colorbar(); c.LineWidth = 1.5; c.FontSize = 20; c.Location = 'southoutside'; c.Color = 'k';
saveas(f2, [pFNS.OPDD2(1 : end-4), '.png']);
