function pPDD = init_entries(pPDD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% CHECKING ENTRIES OF OPDD PARAMETERS
%%%% IF ANY OPTIONAL ENTRIES ARE MISSING, FILL WITH DEFAULT VALUE. WARNING
%%%% WILL RETURN TO NOTIFY THE USERS
%%%% MAIN REFERENCE:
%%%%    [1] LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%% INPUT:
%%%%        PPDD            --> PARAMETERS FOR OPDD 
%%%% OUTPUT:
%%%%        PPDD            --> UPDATED PARAMETERS FOR OPDD 
%%%% CONTACT INFORMATION:
%%%%    STEPHEN Z. LIU 
%%%%    DEPARTMENT OF BIOMEDICAL ENGINEERING
%%%%    JOHNS HOPKINS UNIVERSITY SCHOOL OF MEDICINE
%%%%    E-MAIL: SZLIU29@GMAIL.COM
%%%%    PHONE: +1 657-275-2162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(pPDD, 'UseGPU'); pPDD.UseGPU = true; warning('Setting DEFAULT UseGPU = true'); end
if ~isfield(pPDD, 'UseMomentum'); pPDD.UseMomentum = true; warning('Setting DEFAULT UseMomentum = true'); end
if ~isfield(pPDD, 'UseHuber'); pPDD.UseHuber = false; warning('Setting DEFAULT UseHuber = false'); end

if logical(pPDD.UseMomentum)
    if ~isfield(pPDD, 'MomentMemory'); pPDD.MomentMemory = 100; warning('Setting DEFAULT MomentMemory = 100'); end
    if ~isfield(pPDD, 'MomentTruncate'); pPDD.MomentTruncate = [-inf inf]; warning('Setting DEFAULT MomentTruncate = [-inf inf]'); end
else
    pPDD.MomentMemory = 1; % IN THIS CASE WON'T USE MOMENTUM
end

if logical(pPDD.UseHuber)
    if ~isfield(pPDD, 'HuberImageSize'); error('HuberImageSize is MISSING (i.e., 2D size of radiograph in pixel unit)!!'); end
    if ~isfield(pPDD, 'HuberBeta'); error('HuberBeta is MISSING (i.e., initial regularization strengths)!!'); end
    if ~isfield(pPDD, 'HuberDelta'); error('HuberDelta is MISSING (i.e., huber cutoff)!!'); end
    if ~isfield(pPDD, 'UseHuberAdjust'); pPDD.UseHuberAdjust = true; warning('Setting DEFAULT UseHuberAdjust = true'); end
    if logical(pPDD.UseHuberAdjust)
        if ~isfield(pPDD, 'HuberAdjustRatio'); pPDD.HuberAdjustRatio = 0.4; warning('Setting DEFAULT HuberAdjustRatio = 0.4'); end
        if ~isfield(pPDD, 'HuberAdjustFreq'); pPDD.HuberAdjustFreq = 200; warning('Setting DEFAULT HuberAdjustFreq = 200'); end
        if ~isfield(pPDD, 'HuberAdjustRange'); pPDD.HuberAdjustRange = 10; warning('Setting DEFAULT HuberAdjustRange = 20'); end
    end
end

if ~isfield(pPDD, 'BacktrackInitStep'); pPDD.BacktrackInitStep = 0.1; warning('Setting DEFAULT BacktrackInitStep = 0.1'); end
if ~isfield(pPDD, 'BacktrackWolfeC1'); pPDD.BacktrackWolfeC1 = 0.0001; warning('Setting DEFAULT BacktrackWolfeC1 = 0.0001'); end
if ~isfield(pPDD, 'BacktrackWolfeC2'); pPDD.BacktrackWolfeC2 = 0.9999; warning('Setting DEFAULT BacktrackWolfeC2 = 0.9999'); end
if ~isfield(pPDD, 'DrawConvergence'); pPDD.DrawConvergence = true; warning('Setting DEFAULT DrawConvergence = true'); end
if ~isfield(pPDD, 'DrawUpdateIter'); pPDD.DrawUpdateIter = 100; warning('Setting DEFAULT DrawUpdateIter = 100'); end
if ~isfield(pPDD, 'ConvergeTolerance'); pPDD.ConvergeTolerance = 0.0001; warning('Setting DEFAULT ConvergeTolerance = 0.0001'); end
if ~isfield(pPDD, 'ConvergeMaxIter'); pPDD.ConvergeMaxIter = 2000; warning('Setting DEFAULT ConvergeMaxIter = 2000'); end
if ~isfield(pPDD, 'ConvergeRange'); pPDD.ConvergeRange = 20; warning('Setting DEFAULT ConvergeRange = 20'); end


