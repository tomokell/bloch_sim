% Make an SSFP/SPGR sequence for use in Bloch simulations.
%
% Tom Okell, June 2023
%
%
% [G, RF, Spoil, TEIdx, t, ADC] = make_SSFP_seq(FlipAng, RF_dur, TR, dt, T, G_amp, SeqType, NoLinIncPulses, VFAType, G_s, delta, ADC_dur);
%
% This function creates RF and gradients for a bSSFP/SPGR sequence with the
% given parameters.  RF shape is a Hann-windowed sinc by default.  
%
% Inputs:
%   FlipAng         =   Flip angle (degrees), or an array of flip angles for
%                       VFA schemes (see below)
%   RF_dur          =   Excitation RF pulse duration (s)
%   TR              =   Repetition time (s)
%   dt              =   Time increment for waveform sampling (s)
%   T               =   Total sequence duration (s)
%   G_amp           =   Flat-top gradient amplitude [slice select, readout] (mT/m)
%   SeqType         =   'SPGR': spoiled gradient echo
%                       'NoCat': bSSFP with no catalysation                    
%                       'HalfAngle': bSSFP with the first TR having a half
%                       angle excitation and half the TR to catalyse the
%                       steady state
%                       'LinInc': bSSFP with a series of linearly increasing flip
%                       angles to catalyse the steady state
%   NoLinIncPulses  =   Number of linearly inreasing pulses to use if
%                       SeqType is 'LinInc'
%   VFAType         =   Variable flip angle schedule type:
%                       'CFA': constant flip angle, defined by FlipAng
%                       'Lin': Linearly increasing flip angles from
%                       FlipAng(1) to FlipAng(2)
%                       'Quad': Quadratic increase in flip angle from
%                       FlipAng(1) to FlipAng(2)
%                       'Maintain': Recursive calculation of flip angles
%                       to maintain a constant SPGR signal from a given spin (in
%                       the absence of relaxation), such that the final
%                       flip angle in the readout is given by FlipAng(1).
%   G_s             =   Gradient slew rate (mT/m/ms = T/m/s)
%   delta           =   Excitation slab thickness (m)
%   ADC_dur         =   ADC duration within the TR (s)

function [G, RF, Spoil, TEIdx, t, ADC] = make_SSFP_seq(FlipAng, RF_dur, TR, dt, T, G_amp, ...
                                                       SeqType, NoLinIncPulses, VFAType, ...
                                                       G_s, delta, ADC_dur)

  % Set default parameters
  if nargin < 1; FlipAng        = 40;     end
  if nargin < 2; RF_dur         = 1e-3;   end
  if nargin < 3; TR             = 10e-3;  end
  if nargin < 4; dt             = 0.1e-3; end
  if nargin < 5; T              = 1;      end
  if nargin < 6; G_amp          = [0 0];  end
  if nargin < 7; SeqType        = 'HalfAngle';  end
  if nargin < 8; NoLinIncPulses = 10;     end
  if nargin < 9; VFAType        = 'CFA';  end
  if nargin < 10; G_s           = 200;    end 
  if nargin < 11; delta         = 0.2;    end 
  if nargin < 12; ADC_dur       = TR/3;   end 
      
  % Declare gamma
  g = GetGamma;
    
  % Determine the number of prep pulses
  if strcmp(SeqType,'HalfAngle')
      NPrepPulses = 1;
  elseif strcmp(SeqType,'LinInc')
      NPrepPulses = NoLinIncPulses;
  else
      NPrepPulses = 0;
  end

  % Generate the time array, rounding down to the nearest pulse
  N_pulses = floor( T/TR );
  N_img_pulses = N_pulses - NPrepPulses; % Number of RF pulses use for imaging
  disp([ num2str(N_pulses) ' pulses are being generated for this sequence'])
  t = 0:dt:T;

  % Extract slice-select gradient strength from input parameters
  G_ss = G_amp(1);

  % For the purposes of the sinc calculation, if we want gradients off then
  % set the slice-select gradient to a small nominal value, balanced by
  % having a large slab to excite to give a sensible pulse shape
  if G_ss == 0
      G_ss = 1e-9; delta = 1/G_ss;  
  end
  
  % Create a single RF sinc pulse 
  % NB. use the first flip angle here - we will rescale for VFA later if
  % needed
  phi = 0; % Set the phase of the pulse to be zero
  HannWindow = true; % Hann window the pulse to reduce ringing
  RF_pulse = sinc_RF(FlipAng(1), phi, RF_dur, dt, G_ss, g, delta, HannWindow);
  
  % Check the duration of the RF pulse (some functions add or remove one data
  % point to make the pulse symmetric about its centre
  pulse_samples = size(RF_pulse,2);
  
  % Generate balanced gradient pulses for slice selection on the z axis
  [G_pulse_ss, N_ss, flat_top_start_idx_ss] = balanced_grad_pulses(dt, RF_dur, G_s, G_amp(1), [0 0 1]');
  
  % Generate balanced gradient pulses for readout on the x axis
  [G_pulse_ro, N_ro, flat_top_start_idx_ro] = balanced_grad_pulses(dt, ADC_dur, G_s, G_amp(2), [1 0 0]');

  % Calculate flip angles required for maintaining the signal level by
  % compensating for previous RF pulses
  if strcmp(VFAType,'Maintain')
      MaintainVFAs = MaintainVFA(N_img_pulses,FlipAng(1));
  end
    
  % Initialise outputs
  RF = zeros(2,size(t,2));    G = zeros(3,size(t,2));  Spoil = false(size(t)); ADC = false(size(t));
  TEIdx = zeros(1,N_img_pulses);
  
  % Loop through the TRs, generating the gradient waveforms for each
  TRStartIdx = 1;

  for ii = 1:N_pulses
  
      % Calculate the RF scaling factor required for prep pulses (half
      % angle/linear increase)
      if ii <= NPrepPulses 
          RFcount = 1; % This results in prep pulse scaling relative to the first excitation pulse for VFA sequences
          
          if strcmp(SeqType,'HalfAngle')  
              PrepPulseFactor = 0.5;

          elseif strcmp(SeqType, 'LinInc')
              PrepPulseFactor = ii/(NoLinIncPulses+1);

          else 
              error(['Unknown sequence type: ' SeqType '!'])

          end

      else % Not a prep pulse
          PrepPulseFactor = 1;
          RFcount = ii - NPrepPulses; % Count from the end of the prep pulses now
      end

      % Scale for VFA
      switch VFAType
          case 'CFA' % Constant flip angle
              RFscalefactor = 1;

          case 'Lin' % Linear increase
              RFscalefactor = 1+(RFcount-1)/(N_img_pulses-1)*(FlipAng(2)-FlipAng(1))/FlipAng(1);

          case 'Quad' % Quadratic increase
              RFscalefactor = 1+((RFcount-1)/(N_img_pulses-1))^2*(FlipAng(2)-FlipAng(1))/FlipAng(1);

          case 'Maintain' % Maintains signal level by compensating for previous RF pulses
              RFscalefactor = MaintainVFAs(RFcount)/MaintainVFAs(end);

          otherwise
              error('VFA mode unknown')
      end
          
      % Calculate the RF pulse start/end indices
      pulse_start = flat_top_start_idx_ss + TRStartIdx - 1;
      pulse_end = pulse_start + pulse_samples - 1;

      % Copy the single pulse into the RF array at this index, scaling for
      % VFA if needed
      RF(1, pulse_start:pulse_end ) = RF_pulse(1,:)*RFscalefactor*PrepPulseFactor;
      RF(2, pulse_start:pulse_end ) = RF_pulse(2,:);

      % Determine the start of the next TR
      if (ii == 1) && strcmp(SeqType,'HalfAngle') % Half the TR for the first HalfAngle case
          NextTRStartIdx = TRStartIdx + round(TR/2/dt);
      else
          NextTRStartIdx = TRStartIdx + round(TR/dt);
      end

      % Spoil at the end of the TR if required
      if strcmp(SeqType,'SPGR')
          Spoil(NextTRStartIdx-1) = true;
      end

      % For bSSFP, set the RF phase to alternate between 0 and pi
      if ~strcmp(SeqType,'SPGR')
        RF(2,pulse_start:(pulse_start + pulse_samples - 1) ) = mod(ii+1,2) * pi;
      end

      % Insert the slice-select gradient pulse into the array      
      G(:,TRStartIdx:(TRStartIdx + N_ss - 1) ) = G_pulse_ss;

      % Calculate the readout gradient/ADC timings if we're not in prep mode
      if ii > NPrepPulses
          % Calculate the time available
          EarliestStartIdx = TRStartIdx + flat_top_start_idx_ss + ceil(RF_dur/dt);

          % Determine the TE index (half way between the centre of two RF pulses)
          TEIdx(ii-NPrepPulses) = round((TRStartIdx+NextTRStartIdx)/2 + pulse_samples/2);          
          
          % Put the readout centred on the TE
          ROStartIdx = TEIdx(ii-NPrepPulses)-floor(N_ro/2);
          ROEndIdx = ROStartIdx + N_ro - 1;

          % Check there is sufficient time in the TR
          if ROStartIdx < EarliestStartIdx
              error(['Not enough time within the TR for the readout! Needs ' ns((EarliestStartIdx-ROStartIdx)*2*dt) ' s more!'])
          end

          % Insert the readout gradient pulses (adding to the existing
          % gradients so they can overlap if needed).
          G(:,ROStartIdx:ROEndIdx) = G(:,ROStartIdx:ROEndIdx) + G_pulse_ro;

          % ADC
          ADCStart = ROStartIdx + flat_top_start_idx_ro - 1;
          ADCEnd = ADCStart + ceil(ADC_dur/dt) - 1;
          ADC(ADCStart:ADCEnd) = true;
          
      end

      % Set the start index for the next TR
      TRStartIdx = NextTRStartIdx;
  end
  