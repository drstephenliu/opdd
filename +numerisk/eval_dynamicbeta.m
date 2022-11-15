function HuberBeta = eval_dynamicbeta(iter, HuberBeta, ObjectiveReg, ...
    ObjectiveTotal, ratio, range, frequency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% ADJUST REGULARIZATION STRENGTH BASED ON THE RATIO BETWEEN
%%%% REGULARIZATION AND THE TOTAL OBJECTIVE FUNCTION
%%%% MAIN REFERENCE:
%%%%    [1] LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%% INPUT:
%%%%        ITER            --> CURRENT ITERATION COUNTER
%%%%        HUBERBETA       --> STRENGTH OF PREVIOUS ITERATION
%%%%        OBJECTIVEREG    --> REGULARIZATION VALUE
%%%%        OBJECTIVETOTAL  --> TOTAL OBJECTIVE VALUE
%%%%        RATIO           --> DESIRED RATIO REGULARIZATION/TOTAL
%%%%        RANGE           --> COMPUTED OVER THE MOST RECENT N ITERATION
%%%%        FREQUENCY       --> THE FREQUENCY OF ADJUSTING REGULARIZATION
%%%% OUTPUT:
%%%%        HUBERBETA       --> UPDATED STRENGTH
%%%% CONTACT INFORMATION:
%%%%    STEPHEN Z. LIU 
%%%%    DEPARTMENT OF BIOMEDICAL ENGINEERING
%%%%    JOHNS HOPKINS UNIVERSITY SCHOOL OF MEDICINE
%%%%    E-MAIL: SZLIU29@GMAIL.COM
%%%%    PHONE: +1 657-275-2162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ratio = max(min(ratio, 1), 0);
if frequency < range; error('FREQUENCY OF ADJUSTMENT SMALLER THAN COMPUTED RANGE!'); end

if (rem(iter, frequency) == 0)
    test_obj = mean(ObjectiveReg(iter-(range-1) : iter));
    test_objective = mean(ObjectiveTotal(iter-(range-1) : iter));
    test_scalar = ratio / (test_obj / test_objective);
    
    HuberBeta = HuberBeta * test_scalar;
    warning('BETA AUTO-ADJUSTED!!!');
end