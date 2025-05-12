%-----------------------------------------------------------------------------------------------------
% MATLAB simulation of a basic MUSIC direction finding algorithm
% Mark Patten
% May 12, 2025
% 
% Generates 100 snapshots received by 8 antenna linear array spaced half a meter apart. 
% Simulation is of three transmitters 300 MHz,
%   spaced at 10, 40, and 60 degrees from broadside to the antenna array,
%   each received at 10 dB. 
% 
% The number of transmitters, as well as their azimuth location and power, can easily be   %   changed by modifying the definitions of the “txAngles” and “txSNRs” arrays.
% 
% The code takes the snapshots computed and calculates MUSIC spectrum,
%   showing peaks at 10, 40, and 60 degrees:
%----------------------------------------------------------------------------------------------------

%------------------------------------------
% set up initial simulation parameters
%------------------------------------------
numAnts = 8; % number of antennas
antSpacing = 0.5; % linear array of antennas 0.5 meters apart
frequency = 300e6; % 300 MHz
speedLight = 2.9979e8; % speed of light in meters/second
numSnaps = 100; % number of snapshots
txSNRs = [10 10 10]; % SNRs for three transmitters
txAngles = [10 40 60]; % angles of arrival in degrees for three transmitters
%----------------------------------------------
% generate simulated transmitted signals
%----------------------------------------------
numTxs = length(txAngles); % number of transmitters
wavelen = speedLight/frequency;
txSigs = []; % initialize transmitter signals array to blank
for txIdx = 1:numTxs % loop for each transmitter
   phase = 2*pi*rand([numSnaps 1]); % generate uniform random phase
   txSig = (cos(phase)+i*sin(phase)) ... % generate vector with constant modulus
      *10^(txSNRs(txIdx)/20); % multiply by signal magnitude
   txSigs =[txSigs txSig]; % add this transmitter signal to signals array
end
%---------------------------------------------------
% generate snapshot array of antenna signals
%---------------------------------------------------
antSnaps = []; % initialize snapshot array to blank
for antIdx = 1:numAnts % loop for each antenna
   % initialize snapshot vector with unity variance gaussian noise
   snaps = (sqrt(2)/2)*(randn(numSnaps,1)+i*randn(numSnaps,1));
   for txIdx = 1:numTxs % loop for each transmitter
      % compute phase shift for this transmitter's directon angle
      phaseShift = sind(txAngles(txIdx))*(antIdx-1)*antSpacing/wavelen;
      phaseRotator = cos(phaseShift*2*pi) + i*sin(phaseShift*2*pi);
      % add phase shifted transmitter signal to snaps vector for this antenna
      snaps = snaps + txSigs(:,txIdx)*phaseRotator;
   end
   antSnaps = [antSnaps snaps]; % insert snaps vector into antenna snaps array
end
%------------------------------------------------------
% compute MUSIC direction spectrum
%------------------------------------------------------
covariance=antSnaps'*antSnaps; % compute antenna snapshot covariance
[vx dx]=eig(covariance); % compute eigenvalues and eigenvectors
[eigValues,idx] = sort(diag(dx)); % Vector of sorted eigenvalues
eigVects = vx(:,idx); % Sort eigenvectors from eigenvalues
nsEigVects = eigVects(:,1:end-numTxs); % Noise Space eigenvectors
steeringAngles=-90:1:90;
musicSpec = zeros(size(steeringAngles));
for steeringAngleIdx = 1:length(steeringAngles)
   antPhaseDelays = [0:1:(numAnts-1)]' ...
   *sind(steeringAngles(steeringAngleIdx))*antSpacing/wavelen;
   antSteeringVect = cos(antPhaseDelays*2*pi)-i*sin(antPhaseDelays*2*pi);
   musicSpecVal = 1/(antSteeringVect'*nsEigVects*nsEigVects'*antSteeringVect);
   musicSpec(steeringAngleIdx) = musicSpecVal;
end
%---------------------------
% plot MUSIC spectrum
%---------------------------
plot(steeringAngles,abs(musicSpec))
