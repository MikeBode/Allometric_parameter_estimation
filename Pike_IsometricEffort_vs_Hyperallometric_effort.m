function Pike_Isometric_effort_vs_Hyperallometric_effort()
clc
PT = 1;

%% Calculate the desired steepness, assuming isometry
% Steepness is defined as equilibrium recruitment (expressed as a proportion of virgin recruitment) 
% observed when the stock is reduced to 20% of its virgin biomass.

% Start by assuming that the population has isometric reproduction scaling
ISO = 1;

% Load the parameters for pike
CreateParameters(ISO)

% Set the steepness goal
Steepness_goal = 0.5;

% Choose the log-transformed BH parameters 
% (log transforming makes them numerically better behaved for searches) 
%           R = exp(LN_a)*S / (1 + exp(LB_b)*S))
% where R is the number of recruits entering the first age class,
%       S is the number of settlers, pre-DD mortality.
% These parameters help to set the steepness
% For pike, a steepness of 0.6 results from LN_a = -4.75
%           a steepness of 0.5 results from LN_a = -5.25
%           a steepness of 0.4 results from LN_a = -5.72
%           a steepness of 0.3 results from LN_a = -6.20
LN_a = -5.25;
LN_b = -5; % I'm keeping this constant

% Measure recruitment & biomass across harvest rates to calculate the steepness
H_vec = linspace(0,0.75,250); % These are the fishing rates being considered
UB = 1; % This is the upper bound of the x-axis
Abundance = zeros(size(H_vec)); Yield_biomass = Abundance; Biomass = Abundance; recruits = Abundance;
for h = 1:length(H_vec)
    [Abundance(h),Yield_biomass(h),Biomass(h),recruits(h)] = sub_beverton_holt_model(LN_a,LN_b,H_vec(h),[],ISO);
    if Biomass(h) < 1e-2 % Once the population has crashed, there's no point looking at higher harvest rates
        UB = H_vec(h); % Keep track of this for plotting
        break
    end
end
disp(['Biomass scale = ' num2str(round(Biomass(1)))]) % Display the virgin biomass (e.g., when harvest = 0)
NormBiomass = Biomass./Biomass(1); % Normalise biomass and recruitment to virgin levels
NormRecruits = recruits./recruits(1);
VirginRecruitment_Iso = recruits(1); % Keep track of this to standardise the hyperallometric 

% Verify that we've achieved the target steepness
% Where is biomass 20% of virgin?
[BiomassResidual,I] = min(abs(Biomass - 0.2*Biomass(1)));
Steepness_achieved = recruits(I)./recruits(1); % This is the current steepness
disp(['Steepness = ' num2str(Steepness_achieved,2)])

% Now calculate the harvest mortality that gave MSY (biomass) for the isometric population
[MSY,I_MSY] = max(Yield_biomass);
NormYield_biomass_iso = Yield_biomass./max(Yield_biomass);
F_MSY_Isometric = H_vec(I_MSY);


% Plot the steepness calculation figure
figure(2), clf; 
hold on; FS = 16; CL = get(gca,'colororder');
plot(H_vec,NormBiomass,'-','linewidth',2);
plot(H_vec,NormRecruits,'-','linewidth',2);
plot([H_vec(I) H_vec(I) 0],[0 NormRecruits(I) NormRecruits(I)],'k--')
plot([H_vec(I) 0],[NormBiomass(I) NormBiomass(I)],'k--')

%% Inflate the reproductive output to give the hyperallometric population the same virgin recruitment
ISO = 0; % Switch the model to hyperallometric
CreateParameters(ISO); % Recreate the parameters for a hyperallometric population

% Inflate hyperallometric reproductive output by a (log-transformed) constant
% For pike, a steepness of 0.6 requires a F_h value of -1.100
%           a steepness of 0.5 requires a F_h value of -1.240
%           a steepness of 0.4 requires a F_h value of -1.270
%           a steepness of 0.3 requires a F_h value of -1.295
F_h = -1.240; 

disp('--------------------')
disp(['If we reduce hyperallometric output by ' num2str(100*(1-exp(F_h)),2) '%'])
[~,~,~,VirginRecruitment_Hyp] = sub_beverton_holt_model(LN_a,LN_b,0,F_h,ISO);
disp(['Virgin recruitment (isometric)       = ' num2str(VirginRecruitment_Iso,2)])
disp(['Virgin recruitment (hyperallometric) = ' num2str(VirginRecruitment_Hyp,2)])
disp('--------------------')

%% Calculate MSY for the hyperallometric population
Yield_biomass = zeros(size(H_vec));
for h = 1:length(H_vec)
    [~,Yield_biomass(h)] = sub_beverton_holt_model(LN_a,LN_b,H_vec(h),F_h,ISO);
    if Biomass(h) < 1e-2
        UB = max([UB,H_vec(h)]); % Keep track of this for plotting
        break
    end
end

% Now calculate the harvest mortality that gave MSY (biomass) for the isometric population
[MSY,I] = max(Yield_biomass);
NormYield_biomass_hyp = Yield_biomass./max(Yield_biomass);
F_MSY_Hyperallometric = H_vec(I);
Overfishing_ratio = F_MSY_Hyperallometric/F_MSY_Isometric;
disp(['F_MSY (isometric)       = ' num2str(F_MSY_Isometric,3)])
disp(['F_MSY (hyperallometric) = ' num2str(F_MSY_Hyperallometric,3)])
disp(['Hyperallometric F_MSY is ' num2str(round(100*Overfishing_ratio)) '% of isometric F_MSY'])

% Finish the plotting by showing the two yield curves
plot(H_vec,NormYield_biomass_iso,'-','linewidth',1.5,'color',CL(5,:));
plot(H_vec,NormYield_biomass_hyp,'-','linewidth',1.5,'color',CL(6,:));
plot([F_MSY_Isometric F_MSY_Isometric],[0 1],'--','color',CL(5,:))
plot([F_MSY_Hyperallometric F_MSY_Hyperallometric],[0 1],'--','color',CL(6,:))
xlabel('Fishing mortality','fontsize',FS)
xlim([0 UB+0.05]); ylim([0 1.05])
ylabel('Normalised variables','fontsize',FS)
L = legend('Biomass','Recruitment','','','Yield (iso)','Yield (hyp)'); 
set(L,'fontsize',FS,'box','off','location','eastoutside')
text(UB*0.7,1,['Steepness = ' num2str(Steepness_goal)],'fontsize',FS,'color',CL(7,:))
text(UB*0.7,0.8,['$\frac{F_{MSY}^{hyp}}{F_{MSY}^{iso}}$ = ' num2str(round(100*Overfishing_ratio)) '\%'],'fontsize',FS,'interpreter','latex','color',CL(7,:))

FileName = 'F_MSY_comparison.tiff';
Dimensions = [0 0 20 10];
Resolution = '-r200';
set(gcf, 'paperunits', 'centimeters', 'paperposition', Dimensions)
set(gcf, 'renderer', 'painters')
print('-dtiff',Resolution,FileName)

%% ==========================
%% ==== Population model ==== 
%% ==========================
function [Abundance,Yield_biomass,Biomass,recruits] = sub_beverton_holt_model(LN_a,LN_b,H,F_h,ISO)

% The input variables are:
%    LN_a   - Beverton-Holt parameter A (log transformed for numerical stability of search algorithm)
%    LN_b   - Beverton-Holt parameter B (log transformed for numerical stability of search algorithm)
%    H      - Fishing mortality
%    F_h    - Fecundity modifier for hyperallometry (log transformed for numerical stability of search algorithm)
%    ISO    - Binary flag variable for isometry (ISO = 1) or hyperallometry (ISO = 0)

% The output variables are:
%    Abundance          - The total abundance across the metapopulation
%    Biomass            - The total biomass across the metapopulation
%    Yield_abundance    - The total yield in raw individual numbers
%    Yield_biomass      - The total yield in biomass
%    settlers           - The density of settlers
%    recruits           - The density of post_DD recruits

load ParameterValues mass a A f MSL

% Initialise the population age distribution
n = ones(A,1)./A; 

% Beverton Holt parameters
alpha = exp(LN_a);
beta = exp(LN_b);

% Alter total fecundity for hyperallometric populations
V = 1;
if ISO == 1 V = V*1;
elseif ISO == 0 V = V*exp(F_h);
end

% Set a minimum size limit
H = H.*ones(size(a));
H(1:MSL) = 0;

%% Forward simulate the population to numerical equilibrium
delta = 1; t = 1; N = sum(n);
while delta > 1e-4 & t < 1e6
   
   % The model dynamics follow the sequence defined in Hastings & Botsford (1999 Science)
   % - Reproduction
   % - Natural mortality
   % - Settlement
   % - Harvest
   
   % Larvae are produced and put into a pool (before mortality)
   settlers = V.*sum(f.*n); % Larval Pool with Equal Redistribution
   
   % Recruitment mortality occurs according to a Beverton-Holt model. 
   recruits = alpha.*settlers./(1 + beta.*settlers); 
   
   % Older age classes in all areas experience natural mortality (at rate 1-a)
   n = n.*a;
   
   % Individuals get older
   n(2:end) = n(1:end-1);
   
   % Recruitment occurs
   n(1) = recruits; 
   
   % Harvest occurs & catch is recorded
   Catch = n.*H;
   n = n.*(1-H);
   
   % Moving two-year population record
   N = [N(end) sum(n)];
   
   % Exit the loop if we've reached equilibrium (i.e., the total population has ceased to change)
   if t > 2e3; delta = abs(N(2) - N(1)); end; t = t+1;
end

% Translate the catch vector to biomass
Yield_biomass   = sum(Catch.*mass); % Fishing yield (biomass)

% Calculate total biomass & abundance across the age classes
Biomass = sum(n.*mass);
Abundance = sum(n);

% Return an error message if the population is being constrained by our maximum age
if n(end)./sum(n) > 1e-3 % If the final age class is more than 0.1% of the distribution
    disp('Error: maximum age is too low.')
end

%% ========================
%% ==== Parameter file ==== 
%% ========================

function CreateParameters(ISO)
% Load the parameters for pike

% Define the number of age classes, and a corresponding vector
% Technically this shouldn't be finite, but we'll check later if
% we've constrained it too much. 50 years is normally enough.
A = 50; ages = [1:A];

% These are von Bertalaffy parameters and equation, extracted 
% from the vectors for this particular species
L_inf = 118; % Asymptotic length
K = 0.16; % Rate of increase
a0 = -0.509; % This is a negative constant, indicating that the species are not born with zero length
Length = L_inf.*(1 - exp(-K.*(ages - a0)));

% % Calculate natural mortality and survival based on Charnov et al. 
M_continuous = K*(Length./L_inf).^(-1.5);
a = exp(-M_continuous)'; % Convert from a continuous rate to a discrete time probability
NM = 1-a; % Mortality is the complement of survival

% Calculate mass from length according to the length-weight parameters
mass     = 0.0045.*Length.^3.107; mass = mass';
mass_inf = 0.0045.*L_inf .^3.107; 

% Relative fecundity in each age class. It has a coefficient, but we can subsume it into the free variable of larval mortality
if ISO == 1 z = 1; % Over-ride the exponent with isometry if necessary
elseif ISO == 0 z = 1.15; % Split the difference between the exponents fitted to log-transformed and un-transformed data
end
f = mass.^z; 

% Fecundity is zero before a mimumum length
% For pike, age at maturation for females is length dependent, with age 50% maturation 38 cm, roughly age 2 is first reproduction
First_repro = min(find(Length > 38));
f(1:First_repro-1) = 0;

% Often pike are managed with size limits to prevent catching fish below 50 cm (age 3). 
MSL = 3;

save ParameterValues







