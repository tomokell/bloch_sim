%% Simulate the dependence of PCASL inversion efficiency on velocity
% 
% Tom Okell, May 2023

%% Set sequence default parameters
% Similar to Okell MRM 2010
set_VEPCASL_defaults_gaussianhann;

%% Define other parameters
zmax = 3 / 100; % cm to m conversion
cycles = 1:2; % Tag and control cycles (no vessel-encoding)

%% Define velocities to be tested
vs = (1:200) / 100;

%% Initialise arrays
FinalMz = zeros(length(vs),length(cycles)); IE = FinalMz;

%% Simulate for each velocity and cycle number
% Store the final Mz in each case
for jj = 1:length(cycles)
    cycle = cycles(jj);

    parfor ii=1:length(vs)
        v = vs(ii);

        disp(['Velocity ' ns(ii) ', cycle ' ns(jj)])

        % Run the simulation
        [M, t, P, G, RF] = test_VEPCASL_seq_meanGz(v, zmax, z_offset, meanGz, G_amp, RF_shape, ...
                                                   RF_shape_params, RF_amp, RF_dur, ...
                                                   RF_sep, Pa, Pb, cycle, ...
                                                   Ps, dt, T1,T2, ORFreq);

        % Record the final z magnetisation
        FinalMz(ii,jj) = M(3,end);

        % Calculate the inversion efficiency
        IE(ii,jj) = InvEff(FinalMz(ii,jj),2*zmax/v, zmax/v, T1);

    end
end

%% Plot the inversion efficiency vs. velocity for tag and control
figure
plot(vs*100, IE, 'LineWidth', 2);
xlabel 'Velocity in cm/s'
ylabel 'Inversion Efficiency'
ylim([0 1])
legend({'Tag','Control'})

%% Calculate the equivalent curve assuming laminar flow
% NB. Weighting by velocity gives an average over the volume of the blood
% bolus which passes through the labelling plane, rather than an average
% over the spins present at the labelling plane at any given moment, so is
% more meaningful in terms of the average inversion efficiency.
v_avs = vs(2*vs<max(vs)); % Only calculate for average velocities < 2*v_max to ensure accuracy in the laminar flow averaging
IE_Laminar_v = zeros(length(v_avs), size(IE,2));
IE_Laminar   = zeros(length(v_avs), size(IE,2));
for ii = 1:length(v_avs) % Average velocities
    for jj = 1:size(IE_Laminar_v,2) % Cycles
        v_av = v_avs(ii);
        IE_Laminar_v(ii,jj) = AverageOverLaminarFlowProfileWeightedByV(vs,IE(:,jj),v_av);        
        IE_Laminar  (ii,jj) = AverageOverLaminarFlowProfile(vs,IE(:,jj),v_av);        
    end
end

%% Plot the laminar average, weighted by velocity, along with the original
figure
plot(vs*100, IE, 'linewidth', 2); hold on;
plot(v_avs*100, IE_Laminar_v, 'linewidth', 2); 
xlim([min(v_avs) max(v_avs)]*100)
xlabel 'Velocity in cm/s'
ylabel 'Inversion Efficiency'
title 'Laminar weighting vs. unweighted inversion efficiency'
legend({'Tag','Control','Tag: weighted','Control: weighted'})
ylim([0 1])

%% Calculate the contrast (i.e. overall inversion efficiency, accounting for tag and control)
C = IE(:,1) - IE(:,2);
C_Laminar_v = IE_Laminar_v(:,1) - IE_Laminar_v(:,2);

%% Plot the contrast
figure;
plot(vs*100,C, 'linewidth', 2); hold on;
plot(v_avs*100, C_Laminar_v, 'linewidth', 2);
xlabel 'Velocity in cm/s'
ylabel 'Contrast'
legend('Single velocity','Laminar flow: weighted')
ylim([0 1]); xlim([0 max(v_avs)]*100)
