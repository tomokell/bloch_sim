% This function plots the given sequence in a format for use in
% presentations/papers etc.
%
% Tom Okell, June 2023

% Inputs are: time array, t (s), gradient array, G (3xN, mT/m), RF array,
% RF (2xN [amplitude (mT), phase (rads)]), ADC array, ADC (logical, 1xN),
% an array of time indices to be marked as TE, TEIdx, and an array that
% marks when spoiling is performed, Spoil (logical, 1xN).
%
% plot_seq_for_paper(t, G, RF, ADC, TEIdx, Spoil)

function plot_seq_for_paper(t, G, RF, ADC, TEIdx, Spoil)

  if nargin < 4; ADC   = []; end
  if nargin < 5; TEIdx = []; end
  if nargin < 6; Spoil = []; end

  % Deal with empty parameters
  if isempty(G);   G   = zeros(3,length(t)); end
  if isempty(RF);  RF  = zeros(2,length(t)); end
  if isempty(ADC); ADC = zeros(1,length(t)); end
  
  %% Normalise the data to less than unit range
  Frac = 0.9;
  G = G / range(G(:)) * Frac;
  RF(1,:) = RF(1,:) / range(RF(1,:)) * Frac; % Amplitude
  RF(2,:) = RF(2,:) / range(RF(2,:)) * Frac; % Phase
  ADC = ADC / range(ADC(:)) * Frac;

  %% Add offsets to get plots in different parts of the figure
  ADC     = ADC     + 6;
  RF(1,:) = RF(1,:) + 5;
  RF(2,:) = RF(2,:) + 4;
  G(1,:)  =  G(1,:) + 3;
  G(2,:)  =  G(2,:) + 2;
  G(3,:)  =  G(3,:) + 1;

  %% Plot the data and axis lines
  figure;
  plot(t*1e6, [ADC' RF' G'], 'k', 'LineWidth', 1.5);
  hold on;
  plot([t(1) t(end)]*1e6, [1 1; 2 2; 3 3; 4 4; 5 5; 6 6;]','k')
  ylim([0 7.2]);

  set(gca, 'YTick', 1:7);
  set(gca, 'YTickLabel',{'Gz', 'Gy', 'Gx', 'RF phase','RF', 'ADC', ''})
  set(gca, 'XTick', [])
  set(gca, 'TickLength', [0 0])
  set(gca, 'FontSize', 14)
  
  xlim([min(t) max(t)]*1e6)

  %% Add TEIdx lines if supplied
  if ~isempty(TEIdx)
      for ii = 1:length(TEIdx)
          plot(t(TEIdx(ii))*[1 1]*1e6,[min(ADC) max(ADC)+0.1],'g--', 'LineWidth', 1.5)
      end
  end

  %% Add spoiling lines if supplied
  if ~isempty(Spoil)
      SpoilIdx = find(Spoil);
      for ii = 1:length(SpoilIdx)
          plot(t(SpoilIdx(ii))*[1 1]*1e6,[min(G(:)) max(G(:))+0.1],'r:', 'LineWidth', 1.5)
      end
  end