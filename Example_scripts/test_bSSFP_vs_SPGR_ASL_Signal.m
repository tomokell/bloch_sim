% Run example ASL signal simulations when using a bSSFP vs. spoiled GRE
% readouts.
%
% Tom Okell, 2016, updated June 2023
%
% Reproduces simulation figures from Okell TW et al. NMR in Biomedicine.
% 2016; 29:776-786. https://www.doi.org/10.1002/nbm.3515

%% >>>> Spoiled gradient echo vs. bSSFP (Fig 2) <<<< %%
% Here we simulate the magnetisation that would be produced after ASL
% label and control conditions when subjected to either a spoiled gradient
% echo (SPGR) or balanced steady state free precession (bSSFP) readout to
% demonstrate the signal attenuation that occurs in both cases due to the
% excitation pulses. We examine three regimes: i) small flip angle, long
% TR; ii) short TR, but very small flip angle to match the signal
% attenuation of regime (i); iii) large flip angle, short TR. This shows
% the huge improvement in SNR efficiency that can be obtained with a bSSFP
% readout due to the reduced signal attenuation effect.

% >> Regime (i): small flip angle, long TR
% Spoiled gradient echo (spoiling on)
set_SSFP_defaults; % Set default bSSFP parameters
TR = 12e-3; % TR in s
FlipAng = 8; % Flip angle, degs
SeqType = 'SPGR'; % Sequence type (spoiled gradient echo - assumes perfect spoiling)
RF_dur = 1000e-6; % RF pulse duration (s)
NoLinIncPulses = 0; % Linearly increasing flip angles at the start of the readout (helps catalyse bSSFP readouts, see below)
T = 1.1; % Simulation time (s)

% Run the simulation for both tag and control situations
[C_si, tC_si, Mtag_si, Mcntl_si, ~, ~, ~, t_i] = test_SSFP_seq_ASL_contrast(v, G_amp, FlipAng, T, RF_dur, TR, ...
                                            Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses);

% bSSFP (spoiling off)
SeqType = 'HalfAngle'; % bSSFP with a half angle pulse as the first pulse
FlipAng = 16; % Use twice the flip angle to give (approx) the same transverse magnetisation as SPGR
[C_bi, tC_bi, Mtag_bi, Mcntl_bi, ~, G_bi, RF_bi, ~, TEIdx_bi, ~, ~, ADC_bi] = test_SSFP_seq_ASL_contrast(v, G_amp, FlipAng, T, RF_dur, TR, ...
                                            Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses);

% >> Regime (ii): Equivalent attenuation but more efficient sampling (shorter TR, smaller flip angle)
% SPGR
TR = 4.2e-3; SeqType = 'SPGR'; FlipAng = 5; 
[C_sii, tC_sii, Mtag_sii, Mcntl_sii, ~, ~, ~, t_ii] = test_SSFP_seq_ASL_contrast(v, G_amp, FlipAng, T, RF_dur, TR, ...
                                            Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses);

% bSSFP
SeqType = 'HalfAngle'; FlipAng = 10; 
[C_bii, tC_bii, Mtag_bii, Mcntl_bii] = test_SSFP_seq_ASL_contrast(v, G_amp, FlipAng, T, RF_dur, TR, ...
                                            Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses);
                                        
% Regime (iii): Demanding case - short TR, high flip angle
% SPGR
TR = 4.2e-3; SeqType = 'SPGR'; FlipAng = 20; 
[C_siii, tC_siii, Mtag_siii, Mcntl_siii, ~, ~, ~, t_iii] = test_SSFP_seq_ASL_contrast(v, G_amp, FlipAng, T, RF_dur, TR, ...
                                            Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses);
                                        
% bSSFP
SeqType = 'HalfAngle'; FlipAng = 40; 
[C_biii, tC_biii, Mtag_biii, Mcntl_biii] = test_SSFP_seq_ASL_contrast(v, G_amp, FlipAng, T, RF_dur, TR, ...
                                            Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses);
                                        
% Repeat for zero RF power as a reference with no signal attenuation
[C_0, tC_0, Mtag_0, Mcntl_0, ~, ~, ~, t_0] = test_SSFP_seq_ASL_contrast(v, G_amp, 0, T, RF_dur, TR, ...
                                            Ps, dt, T1, T2, ORFreq, SeqType, NoLinIncPulses);
                                        
%% Plot part of the regime 1 bSSFP pulse sequence
TR = 12e-3; Idx = round(TR*5/dt); % Show the first five TRs
plot_seq_for_paper(t_i(1:Idx), G_bi(:,1:Idx), RF_bi(:,1:Idx), ADC_bi(1:Idx), TEIdx_bi(TEIdx_bi<Idx))
title 'bSSFP pulse sequence: regime (i)'
% Note that the green dashed lines indicate the TE and phase encoding
% gradients haven't been implemented in this simulation

%% Plot
LargeFigWindow(0.8,0.6); % Set up figure window

% >> Regime (i): Mz
DeltTIdx = round(12e-3/dt); % Plot the signal once per TR
subplot(2,3,1); 
tIdx = 1:DeltTIdx:length(t_i); % Sample points once per TR
plot(t_i(tIdx),Mtag_si(3,tIdx), 'r--','LineWidth',2); hold on % SPGR tag Mz
plot(t_i(tIdx),Mcntl_si(3,tIdx),'r-','LineWidth',2); % SPGR control Mz
plot(t_i(tIdx),Mtag_bi(3,tIdx) ,'b--','LineWidth',2); % bSSFP tag Mz
plot(t_i(tIdx),Mcntl_bi(3,tIdx),'b-','LineWidth',2);  % bSSFP control Mz
plot(t_0(tIdx),Mtag_0(3,tIdx)  ,'k--'); % No RF tag
plot(t_0(tIdx),Mcntl_0(3,tIdx) ,'k-');  % No RF control
grid on; ylabel 'M_z'; ylim([-1.1 1.1]); xlim([0 1]); 
title('Regime (i): \alpha = 8^{\circ}/16^{\circ}, TR = 12 ms')

legend('SPGR: Tag','SPGR: Control','bSSFP: Tag','bSSFP: Control',...
       'No RF: Tag','No RF: Control','Location',[0.1533    0.7421    0.1084    0.1058])

% >> Regime (ii): Mz
subplot(2,3,2);
DeltTIdx = round(4.2e-3/dt); % Sampling index for shorter TR
tIdx = 1:DeltTIdx:length(t_ii); % Sample points once per TR
plot(t_ii(tIdx),Mtag_sii(3,tIdx), 'r--','LineWidth',2); hold on % SPGR tag Mz
plot(t_ii(tIdx),Mcntl_sii(3,tIdx),'r-','LineWidth',2); % SPGR control Mz
plot(t_ii(tIdx),Mtag_bii(3,tIdx) ,'b--','LineWidth',2); % bSSFP tag Mz
plot(t_ii(tIdx),Mcntl_bii(3,tIdx),'b-','LineWidth',2); % bSSFP control Mz
plot(t_0(tIdx),Mtag_0(3,tIdx)  ,'k--'); % No RF tag
plot(t_0(tIdx),Mcntl_0(3,tIdx) ,'k-'); % No RF control
grid on; ylim([-1.1 1.1]); xlim([0 1]); 
title('Regime (ii): \alpha = 5^{\circ}/10^{\circ}, TR = 4.2 ms')

% >> Regime (iii): Mz
subplot(2,3,3);
DeltTIdx = round(4.2e-3/dt); tIdx = 1:DeltTIdx:length(t_i); % as above
plot(t_iii(tIdx),Mtag_siii(3,tIdx), 'r--','LineWidth',2); hold on % SPGR tag Mz
plot(t_iii(tIdx),Mcntl_siii(3,tIdx),'r-','LineWidth',2); % SPGR control Mz
plot(t_iii(tIdx),Mtag_biii(3,tIdx) ,'b--','LineWidth',2); % bSSFP tag Mz
plot(t_iii(tIdx),Mcntl_biii(3,tIdx),'b-','LineWidth',2); % bSSFP control Mz
plot(t_0(tIdx),Mtag_0(3,tIdx)  ,'k--'); % No RF tag
plot(t_0(tIdx),Mcntl_0(3,tIdx) ,'k-'); % No RF control 
grid on; ylim([-1.1 1.1]); xlim([0 1]); 
%title 'Flip Angle = 20^{\circ}, TR = 4.2 ms'
title('Regime (iii): \alpha = 20^{\circ}/40^{\circ}, TR = 4.2 ms')

% >> Regime (i): ASL signal
subplot(2,3,4)
plot(tC_si,C_si,'r-','LineWidth',2); hold on % SPGR
plot(tC_bi,C_bi,'b-','LineWidth',2); xlim([0 1]); % bSSFP
grid on; xlabel 'Time/s'; ylabel 'ASL Signal'; ylim([0 0.75])
legend('SPGR','bSSFP','Location','NorthWest')

% >> Regime (ii): ASL signal
subplot(2,3,5)
plot(tC_sii,C_sii,'r-','LineWidth',2); hold on % SPGR
plot(tC_bii,C_bii,'b-','LineWidth',2); xlim([0 1]); ylim([0 0.75]) % bSSFP
grid on; xlabel 'Time/s'; 

% >> Regime (iii): ASL signal
subplot(2,3,6)
plot(tC_siii,C_siii,'r-','LineWidth',2); hold on % SPGR
plot(tC_biii,C_biii,'b-','LineWidth',2); xlim([0 1]); ylim([0 0.75]) % bSSFP
grid on; xlabel 'Time/s'; 

% Reposition subplots to save space
for row = 1:2
    for col = 1:3
        % Select the subplot
        subplotno = (row-1)*3 + col;
        subplot(2,3,subplotno)
        
        % Reduce the space in between plots
        pos = get(gca,'Position');%[left bottom width height]
        pos(1) = pos(1) - (col-1)*0.03;
        pos(2) = pos(2) + (row-1)*0.07;
        set(gca,'Position',pos);
    end
end

%% >>>>>>>>> Transient oscillation suppression (Fig. 3) <<<<<<<<< %%
% Here we examine the sensitivity of bSSFP readouts to signal oscillations
% when the magnetisation is not on-resonance. This can partially suppressed
% through the use of a series of linearly increasing flip angles at the
% start of the readout, compared to a conventional half angle initial
% excitation pulse.
clear;

% Set default parameters
set_SSFP_defaults; TR = 4.2e-3; Alpha = 40; RF_dur = 1200e-6; T = 1.1;

%% Compare half-angle and linearly increasing flip angle catalysation schemes at different off-resonance frequencies
ORFreqs = [0 25 50 75 100]; % Off-resonance frequencies to test (Hz)
NoLinIncPulses = [10 20 40]; % Number of linearly increasing flip angles to test
tC_halfang = {}; C_halfang = {}; tC_lininc = {}; C_lininc = {}; % Initialise outputs
Mtag_halfang = {}; Mcntl_halfang = {}; Mtag_lininc = {}; Mcntl_lininc = {};

% Loop through off-resonance frequencies
for ii = 1:length(ORFreqs)
   ORFreq = ORFreqs(ii);
   % Simulate the ASL signal (contrast) for the half angle case
   [C_halfang{ii}, tC_halfang{ii}, Mtag_halfang{ii}, Mcntl_halfang{ii}] = ...
       test_SSFP_seq_ASL_contrast(v, G_amp, Alpha, T, RF_dur, TR, Ps, dt, T1, T2, ORFreq, 'HalfAngle', []);

   % Loop through the protocols with different numbers of linearly
   % increasing pulses
   for jj = 1:length(NoLinIncPulses)
       % Simulate the ASL contrast
       [C_lininc{ii,jj}, tC_lininc{ii,jj}, Mtag_lininc{ii,jj}, Mcntl_lininc{ii,jj}] = ...
           test_SSFP_seq_ASL_contrast(v, G_amp, Alpha, T, RF_dur, TR, Ps, dt, T1, T2, ORFreq, 'LinInc', NoLinIncPulses(jj));
   end
end

%% Plot
LargeFigWindow(0.25,0.8)

% Plot the half angle ASL signals across all off-resonance frequencies
% NB. plot in reverse order to make it easier to see the zero frequency
% curve
subplot(1,1+length(NoLinIncPulses),1); 
for ii = length(C_halfang):-1:1
    plot(tC_halfang{ii},C_halfang{ii}); hold on;
end
title('Half angle','FontWeight','normal')
xlim([0 1]); ylim([0 1.5]); grid on; 
xlabel 'Time/s'; ylabel 'ASL Signal'; 

% Now repeat for varying numbers of linearly increasing flip angle pulses
for jj = 1:length(NoLinIncPulses)
    subplot(1,1+length(NoLinIncPulses),jj+1);

    for ii = length(C_halfang):-1:1
        plot(tC_lininc{ii,jj},C_lininc{ii,jj}); hold on;
    end
    xlim([0 1]); ylim([0 1.5]); xlabel 'Time/s'; grid on; 
    title(['N_p = ' ns(NoLinIncPulses(jj))],'FontWeight','normal')
end

% Plot the off-resonance frequencies as the legend (also in reverse order)
legend([ns(flip(ORFreqs)') repmat(' Hz',length(ORFreqs),1)]); 






%% >>>>>> Compare different bSSFP flip angles (Fig. 4) <<<<<< %%
% Here we simulate a set of different flip angles with a bSSFP readout to
% balance initial signal strength against rate of signal attenuation
clear;

% Set parameters
set_SSFP_defaults; TR = 4.2e-3; RF_dur = 1200e-6; T = 1.1;
Alphas = [10 20 40 60 80]; % Flip angles to simulate
NoLinIncPulses = 20; % Use 20 linearly increasing pulses to catalyse the steady-state
ORFreq = 0; % Assume on-resonance for this simulation
tC_lininc = []; C_lininc = []; % Initialise outputs

% Loop through flip angles to test
for ii = 1:length(Alphas)
   Alpha = Alphas(ii);
   [C_lininc(:,ii), tC_lininc] = test_SSFP_seq_ASL_contrast(v, G_amp, Alpha, T, RF_dur, TR, Ps, dt, T1, T2, ORFreq, 'LinInc', NoLinIncPulses);
end

%% Plot
figure; 
plot(tC_lininc,C_lininc);   
legend([ns(Alphas') repmat('^\circ',length(Alphas),1)]); 
xlim([0 1]); ylim([0 1.5]); xlabel 'Time/s'; grid on;
ylabel 'ASL Signal'



%% >>>>>> Bonus: simulate the bSSFP slice profile <<<<<< %%
% Here we simulate the excitation profile of the bSSFP pulse sequence
% used above, to sanity check the slab width and flip angle are correct.
set_SSFP_defaults; % Set default bSSFP parameters
v = 0; % Simulate zero velocity spins for this
SeqType = 'NoCat'; % Use a non-catalysed sequence so the first pulse is a standard pulse
FlipAng = 40; % simulate a 40 degree flip angle
%dt = 1e-6; % You can uncomment and change dt from its default value here
             % to see the effect on the simulated profile

% Simulate one TR only
T = TR;

% Simulate a range of spin positions
zs = -0.3:0.002:0.3;

% Initialise output Mz
Mz = zeros(size(zs));

% Simulate each spin position
for ii = 1:length(zs)
    [M, t, P, G, RF, TEIdx] = test_SSFP_seq(v, Ms, G_amp, FlipAng, T, RF_dur, TR, ...
	                  [0 0 zs(ii)], dt, T1, T2, ORFreq, SeqType, NoLinIncPulses);

    % Save the Mz value at TE
    Mz(ii) = M(3,TEIdx(1));
end

% Convert to effective flip angle
EffFlipAng = torad2deg(acos(Mz));

% Plot
figure; plot(zs,EffFlipAng,'Linewidth',2); 
xlabel 'z/m'; ylabel 'Effective flip angle/\circ'
title 'bSSFP slice profile'

% NB. This gives the expected flip angle in the centre of the profile and
% the correct width (at the moment the default of 20 cm is used, but this
% can be passed as an argument to make_bSSFP_seq.m). Also note there is
% some dependence on dt - using a smaller dt better samples the pulse
% profile giving a more accurate result, but increases simulation time.