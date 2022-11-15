function varargout = init_momentum(pie, DoGPU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% ALLOCATE SPACE FOR MOMENTUM PARAMETERS
%%%% MAIN REFERENCE:
%%%%    [1] LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%% INPUT:
%%%%        PIE             --> VARIABLE 
%%%%        DOGPU           --> IF USING GPU (LOGICAL/BINARY)
%%%% CONTACT INFORMATION:
%%%%    STEPHEN Z. LIU 
%%%%    DEPARTMENT OF BIOMEDICAL ENGINEERING
%%%%    JOHNS HOPKINS UNIVERSITY SCHOOL OF MEDICINE
%%%%    E-MAIL: SZLIU29@GMAIL.COM
%%%%    PHONE: +1 657-275-2162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout{1} = pie;
if logical(DoGPU)
    varargout{2} = gpuArray(0);
    varargout{3} = gpuArray(0);
    varargout{4} = gpuArray(zeros(size(pie)));
else
    varargout{2} = 0;
    varargout{3} = 0;
    varargout{4} = zeros(size(pie));
end

