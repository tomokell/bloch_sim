%% Simulate the vessel-encoded PCASL spatial modulation function across velocities
% 
%  Tom Okell, May 2023
%
%  First we use parameters to match Wong MRM 2007 (DOI 10.1002/mrm.21293)
%  Then Okell MRM 2010 (DOI: 10.1002/mrm.22458)
%  Then unipolar VEPCASL from Wong and Guo, MAGMA 2012 (DOI: 10.1007/s10334-011-0302-7)

%% Declare the gyromagnetic ratio, gamma
g = GetGamma;  % rad/s/mT            

%% Define velocities and positions
v = [5 10 20 40] / 100; %  cm/s to m/s conversion
zmax = 3 / 100; % cm to m conversion (maximum spin offset from the labelling plane to simulate)
Pa = [-0.05 0]'; Pb = [0 0]'; % Position of vessels "A" and "B" within the labelling plane in m

%% Define the off resonance frequency
ORFreq = 0; % No off-resonance in this simulation

%% Define spin positions to be tested
Ps = -0.05:0.001:0.05; Ps(2,:) = 0; 
% NB. This moves the test position from vessel "A" to vessel "B" and
% continues beyond it to see one full period of the inversion profile.

%% Define cycle number and z_offset
cycle = 3; % In this cycle we label at the position corresponding to vessel A and control at vessel B
z_offset = 0; % Position the labelling at isocentre for simplicity

%% Define the gradient parameters
meanGz = 0.08 * 1e-4 * 100 * 1000; % Gauss per cm to mT per m conversion:
                                   % Wong = 0.08 G/cm
G_amp = 0.6 * 1e-4 * 100 * 1000;   % Gauss per cm to mT per m conversion
                                   % Wong = 0.6 G/cm

%% Define T1 and T2
T1 = 1.65; % Wong MRM 2007 assumed infinite T1, but doesn't matter too much as we will correct below
T2 = 0.2;  

%% Define RF parameters for the PCASL sequence
RF_amp = 0.04 * 1e-4 * 1000; % Gauss to mT conversion
                             % Wong = 0.04 G
RF_shape = 'hanning';  
RF_shape_params = []; % No extra parameters needed for Hanning pulses
RF_dur = 600e-6;  % RF duration
RF_sep = 960e-6;  % RF "separation" (interval)

%% Define dt to sample RF adequately
dt = RF_dur/20;

%% For each velocity and spin position, calculate the final Mz
FinalMz = zeros(size(Ps,2),size(v,2)); IE = FinalMz;
for jj = 1:size(v,2)
  parfor ii = 1:size(Ps,2)
    disp(['Velocity step ' num2str(jj) ', position step ' num2str(ii)])
    
    % Run the simulation
    [M, t, P, G, RF] = test_VEPCASL_seq_meanGz(v(jj), zmax, z_offset, meanGz, G_amp, RF_shape, ...
                                                RF_shape_params, RF_amp, RF_dur, ...
                                                RF_sep, Pa, Pb, cycle, ...
                                                Ps(:,ii), dt, T1,T2, ORFreq);

    % Record the final z magnetisation
    FinalMz(ii,jj) = M(3,end);

  end    

  % Calculate the Inversion Efficiency
  IE(:,jj) = InvEff(FinalMz(:,jj),2*zmax/v(jj), zmax/v(jj), T1);

end


%% Plot the results
% Phase accrual between two PCASL pulses:
ph_alt = (Ps(1,:) - Pa(1)) / (Pb(1) - Pa(1)) * pi;
figure
plot(ph_alt, IE, 'LineWidth', 2);
xlabel 'Phase alternation/radians'
ylabel 'Inversion Efficiency'
xlim([0 1]*2*pi)
legend( num2str(100*v') ,'location','best' ); % Legend is velocities in cm/s

%% Equivalent plot to Wong MRM 2007 Fig2b:
figure
plot(ph_alt, (0.5-IE)*2, 'LineWidth', 2);
xlim([0 1]*2*pi)
xlabel 'Phase alternation/radians'
ylabel 'Mz without T1 decay'
legend( num2str(100*v') ,'location','best');




%% ---- Now use parameters similar to Okell MRM 2020 ----- %%
% Set default parameter values
set_VEPCASL_defaults_gaussianhann;

%% Define the velocities to test
% NB. Use a larger range here so we can perform laminar flow averaging
v = [1 3 5 10 20 30 40 60 80 100] / 100; %  cm/s to m/s conversion

%% Set up maximum z value to simulate
zmax = 0.03;

%% Define relative spin positions to be tested
PsRel = -1:0.01:1; 

%% For each velocity and spin position, calculate the final Mz
FinalMz = zeros(size(PsRel,2),size(v,2)); IE = FinalMz;

for jj = 1:size(v,2)
  parfor ii = 1:size(PsRel,2) % Use parfor for speed here, but can be replaced by standard "for" if needed
    disp(['Velocity step ' ns(jj) ' of ' ns(size(v,2)) ', position step ' ns(ii) ' of ' ns(size(PsRel,2))])
    
    % Calculate the absolute spin position
    Ps = Pb + PsRel(ii)*(Pa-Pb);

    % Run the simulation
    [M, t, P, G, RF, T] = test_VEPCASL_seq_meanGz(v(jj), zmax, z_offset, meanGz, G_amp, RF_shape, ...
                                           RF_shape_params, RF_amp, RF_dur, ...
                                           RF_sep, Pa, Pb, cycle, Ps, dt, T1, ...
                                           T2, ORFreq, Unipolar);


    % Record the final z magnetisation
    FinalMz(ii,jj) = M(3,end);

    % Calculate the Inversion Efficiency
    IE(ii,jj) = InvEff(FinalMz(ii,jj),T, T/2, T1);       
    
  end
end


%% Calculate the average over a laminar flow profile, with and without weighting by velocity
% NB. Weighting by velocity gives an average over the volume of the blood
% bolus which passes through the labelling plane, rather than an average
% over the spins present at the labelling plane at any given moment, so is
% more meaningful in terms of the average inversion efficiency.
v_avs = v(2*v<max(v)); % Only calculate for average velocities < 2*v_max to ensure accuracy in the laminar flow averaging
IE_Laminar_v = zeros(size(IE,1),length(v_avs));
IE_Laminar   = zeros(size(IE,1),length(v_avs));
for ii = 1:size(IE_Laminar_v,1)
    for jj = 1:size(IE_Laminar_v,2)
        v_av = v_avs(jj);
        IE_Laminar_v(ii,jj) = AverageOverLaminarFlowProfileWeightedByV(v,IE(ii,:),v_av);        
        IE_Laminar  (ii,jj) = AverageOverLaminarFlowProfile(v,IE(ii,:),v_av);        
    end
end

%% Plot single velocity inversion efficiency
ph_alt = (PsRel + 1) * pi;
figure
plot(ph_alt, IE, 'LineWidth', 2);
xlabel 'Phase alternation/radians'
ylabel 'Inversion Efficiency'
title 'Single velocity (plug flow)'
xlim([0 1]*2*pi)
legend( num2str(100*v'),'location','best' ); % Legend is velocities in cm/s

%% Plot laminar flow averaged version, weighting by v
figure
plot(ph_alt, IE_Laminar_v, 'LineWidth', 2);
xlabel 'Phase alternation/radians'
ylabel 'Inversion Efficiency'
title 'Laminar flow, weighted by v'
xlim([0 1]*2*pi)
legend( num2str(100*v_avs'),'location','best' ); % Legend is velocities in cm/s

%% Equivalent plot to Okell MRM 2010 Fig. A1
figure
psi = mod(ph_alt - pi/2 + pi, 2*pi) - pi; % Recalculate relative distance from encoding centre, as per Chappell MRM 2010
[psi_s,Idx] = sort(psi); % Resort based on this index
plot(psi_s, (0.5-IE_Laminar(Idx,v_avs == 0.3))*2, 'k-','LineWidth', 2); % Plot effective Mz without T1 decay
xlabel 'Psi'
ylabel 'Mz without T1 decay'
xlim([-1 1]*pi); ylim([-1 1]*1.1)
% NB. Standard laminar flow averaging was used here, rather than
% weighting by velocity, which was probably suboptimal.



%% ---- Repeat for unipolar VEPCASL ----- %%
% Ref: Wong and Guo, MAGMA 2012, DOI: 10.1007/s10334-011-0302-7
% Set parameters
g = GetGamma;  % rad/s/mT            
vs = (5:5:40) / 100; %  cm/s to m/s conversion
zmax = 2 / 100; % cm to m conversion
Pa = [0 0]'; Pb = [0.1 0]'; % Position of vessels "A" and "B" within the labelling plane
% NB. Here we put the tag location at the centre, to recreate Fig 1 of the
% paper above.
Ps = -Pb(1):(2*Pb(1)/100):Pb(1); Ps(2,:) = 0; % Evenly space spin test positions from -Pb to Pb
ORFreqs = [0 200]; % Simulate on-resonance and 200 Hz off-resonance
cycles = 3:4; % In these cycles we label at the position corresponding to vessel A and control at vessel B and vice versa
z_offset = 0; % Position the labelling at isocentre for simplicity
meanGz = 0.8; % mT per m. NB. This differs from the paper above, but gives a better match
G_amp = 6; % mT per m
T1 = 100; % Assume a large value
T2 = 0.2;  
RF_shape = 'hamming';  
RF_shape_params = []; % No extra parameters needed for Hamming pulses
RF_dur = 800e-6;  % RF duration
RF_sep = 1400e-6;  % RF "separation" (interval)
dt = RF_dur/40; % Define dt to be smaller here to better represent the effect of off-resonance during RF pulses
RF_amp = 0.04 * 1e-4 * 1000; % Assume the same as Wong MRM 2007
% NB. We can rescale the RF amplitude to achieve a particular flip angle
% like this if we wish:
% DesiredFlipAngle = 20; % Rescale to achieve the desired flip angle
% RF_amp = RF_amp * DesiredFlipAngle / FlipAngle(Hamming_window(RF_dur/dt)*RF_amp,dt);
Unipolar = true; % Unipolar VEPCASL is on (default is bipolar)

%% For each off-resonance, encoding cycle, velocity and spin position, calculate the final Mz
FinalMz = zeros(size(Ps,2),length(vs),length(cycles),length(ORFreqs)); IE = FinalMz;

for nn = 1:length(ORFreqs)
    for kk = 1:length(cycles)
        for jj = 1:length(vs)
            parfor ii = 1:size(Ps,2)
                disp(['Off-resonance step ' num2str(nn) ', cycle ' num2str(kk) ', velocity step ' ns(jj) ', position step ' num2str(ii)])

                % Run the simulation
                [M, t, P, G, RF] = test_VEPCASL_seq_meanGz(vs(jj), zmax, z_offset, meanGz, G_amp, RF_shape, ...
                    RF_shape_params, RF_amp, RF_dur, ...
                    RF_sep, Pa, Pb, cycles(kk), ...
                    Ps(:,ii), dt, T1,T2, ORFreqs(nn), Unipolar);

                % Record the final z magnetisation
                FinalMz(ii,jj,kk,nn) = M(3,end);

            end

            % Calculate the Inversion Efficiency
            IE(:,jj,kk,nn) = InvEff(FinalMz(:,jj,kk,nn),2*zmax/vs(jj), zmax/vs(jj), T1);

        end
    end
end

%% Plot the results for a single velocity
% Phase accrual between two PCASL pulses:
ph_alt = (Ps(1,:) - Pa(1)) / (Pb(1) - Pa(1)) * pi;
ExampleVelocity = 20/100;
figure; hold on; leg = {};
for nn = 1:length(ORFreqs)
    for kk = 1:length(cycles)
        for jj = find(vs==ExampleVelocity)
            plot(ph_alt, IE(:,jj,kk,nn), 'LineWidth', 2);
            leg{end+1} = ['Off-resonance = ' ns(ORFreqs(nn)) ' Hz, cycle = ' ns(cycles(kk)) ', velocity = ' ns(vs(jj)*100) ' cm/s'];
        end
    end
end
xlabel 'Phase alternation/radians'
ylabel 'Inversion Efficiency'
xlim([-1 1]*pi)
legend( leg ,'location','northoutside' ); % Legend is off-resonance/Hz

%% Average over a laminar flow profile
% NB. Weighting by velocity gives an average over the volume of the blood
% bolus which passes through the labelling plane, rather than an average
% over the spins present at the labelling plane at any given moment, so is
% more meaningful in terms of the average inversion efficiency.
v_avs = 20/100; % Equivalent to averaging over velocities up to 40 cm/s
IE_Laminar_v = zeros(size(IE,1),length(v_avs),size(IE,3),size(IE,4));
IE_Laminar   = zeros(size(IE,1),length(v_avs),size(IE,3),size(IE,4));

for ii = 1:size(IE,1)
    for nn = 1:length(ORFreqs)
        for kk = 1:length(cycles)
            for jj = 1:length(v_avs)
                v_av = v_avs(jj);
                IE_Laminar_v(ii,jj,kk,nn) = AverageOverLaminarFlowProfileWeightedByV(vs,IE(ii,:,kk,nn),v_av);
                IE_Laminar  (ii,jj,kk,nn) = AverageOverLaminarFlowProfile(vs,IE(ii,:,kk,nn),v_av);
            end
        end
    end
end

%% Equivalent plot to Wong and Guo MAGMA 2012 Fig 1:
Contrast_Laminar = squeeze((IE_Laminar(:,:,1,:) - IE_Laminar(:,:,2,:)))*2;
Contrast_Laminar_v = squeeze((IE_Laminar_v(:,:,1,:) - IE_Laminar_v(:,:,2,:)))*2;
figure
plot(ph_alt, Contrast_Laminar, 'LineWidth', 2);
xlim([-1 1]*pi)
xlabel 'Phase alternation/radians'
ylabel 'Tagging efficiency (Control-Tag)'
legend( num2str(ORFreqs') ,'location','best');
% NB. This is similar but not identical to Wong and Guo, likely due to some
% of the assumptions made about RF amplitude etc. above.
