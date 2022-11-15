function [pie, pie0, time, time_sum, wind] = eval_momentum(iter, deltas, ...
    pie, pie0, time, time_sum, wind, memorization, truncate_low, truncate_high)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2022/11/10
%%%% NESTEROV MOMENTUM UPDATE GIVEN THE DESCENT STEP
%%%% MAIN REFERENCE:
%%%%    [1] LIU ET AL (2022). "MODEL-BASED THREE-MATERIAL DECOMPOSITION IN
%%%%        DUAL-ENERGY CT USING THE VOLUME CONSERVATION CONSTRAINT". PHYS.
%%%%        MED. BIOL., 64, 145006. DOI: 10.1088/1361-6560/AC7A8B
%%%%    [2] KIM ET AL (2014). "COMBINING ORDERED SUBSETS AND MOMENTUM FOR
%%%%        ACCELERATED X-RAY CT IMAGE RECONSTRUCTION". IEEE TRANS. MED.
%%%%        IMAG., 34(1), 167-178. DOI: 10.1109/TMI.2014.2350962
%%%% INPUT:
%%%%        ITER            --> CURRENT ITERATION COUNTER
%%%%        DELTAS          --> UPDATE DESCENT STEP
%%%%        PIE             --> VARIABLE OF PREVIOUS ITERATION
%%%%        PIE0            --> VARIABLE OF INITIAL MOMENTUM
%%%%        MEMORIZATION    --> LENGTH OF MOMENTUM BEFORE RESET
%%%%        TRUNCATE_LOW    --> LOWER BOUND OF VARIABLE (-INF IF NOT USED)
%%%%        TRUNCATE_HIGH   --> UPPER BOUND OF VARIABLE (+INF IF NOT USED)
%%%%        TIME, TIME_SUM, WIND (PARAMETERS OF PREVIOUS ITERATION)
%%%% OUTPUT:
%%%%        PIE             --> VARIABLE OF CURRENT ITERATION
%%%%        PIE0            --> VARIABLE OF UPDATED INITIAL MOMENTUM
%%%%        TIME, TIME_SUM, WIND (UPDATED PARAMETERS OF CURRENT ITERATION)
%%%% CONTACT INFORMATION:
%%%%    STEPHEN Z. LIU 
%%%%    DEPARTMENT OF BIOMEDICAL ENGINEERING
%%%%    JOHNS HOPKINS UNIVERSITY SCHOOL OF MEDICINE
%%%%    E-MAIL: SZLIU29@GMAIL.COM
%%%%    PHONE: +1 657-275-2162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

memorization = max(memorization, 1);

if memorization == 1 % CONVENTIONAL UPDATE IF MEMORIZATION IS 1 OR LESS
    pie = pie - deltas; 
else
    if rem(iter, memorization) == 0 % RE-INITIALIZE MOMENTUM IF MEMORIZATION IS REACHED
        time = time * 0.0;
        time_sum = time_sum * 0.0;
        wind = wind * 0.0;
        pie0 = pie;
    end
    
    time = 1/2 * (1 + sqrt(1 + 4 * time.^2));
    time_sum = time_sum + time;
    zeta = pie - deltas;
    wind = wind + time * deltas;
    velocity = pie0 - wind;
    pie = zeta + (time/time_sum) * (velocity - zeta);
end

pie = max(min(pie, truncate_high), truncate_low); % TRUNCATE IF WANTED