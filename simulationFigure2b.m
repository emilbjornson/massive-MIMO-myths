%This Matlab script can be used to generate Figure 2b, in the article:
%
%Emil Björnson, Erik G. Larsson, Thomas L. Marzetta, "Massive MIMO: Ten
%Myths and One Critical Question," IEEE Communications Magazine, vol. 54, 
%no. 2, pp. 114-123, February 2016. 
%
%Download article: http://arxiv.org/pdf/1503.06854
%
%This is version 1.01 (Last edited: 2016-04-27)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Initialization
close all;
clear all;


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));

%Select the number of realizations in the Monte-Carlo simulations
nbrOfRealizations = 100000;

%Define the signal-to-noise ratio (SNR)
SNR = 10^(-5/10); % -5 dB

%Define system dimensions
Mmax = 100; %Maximal number of service antennas
K = 12; %Number of terminals
Lmax = 50; %Maximal codebook size


%%Simulation for array gain with Massive MIMO uplink signaling

%Compute estimation quality according to the formula provided in the
%paragraph after Eq. (1)
estimationQuality = 1/(1+1/(K*SNR));

%Compute the array gain as the product between the number of antennas and
%the estimation quality
arrayGainMassiveMIMO = estimationQuality*(1:Mmax);



%%Monte-Carlo simulations for open-loop beamforming with LoS propagation

%Generate uniformly distributed random user terminal angles for LoS
%propagation; note that -pi/2 to pi/2 is analytically equivalent to -pi to
%pi since a ULA cannot separate all angles
userAngles = pi*rand(1,nbrOfRealizations)-pi/2;

%Define distance between antennas in the ULA, measured in wavelengths
interAntennaDistance = 0.5;

%Generate channel realizations with LoS propagation
H_LoS =  exp( repmat((0:Mmax-1)',[1 nbrOfRealizations]) .* repmat(-2i*pi*interAntennaDistance*sin(userAngles),[Mmax 1]) );

%Prepare to compute different array gains
arrayGainLoS = zeros(size(arrayGainMassiveMIMO));


%Go through the range of number of service antennas
for m = 1:Mmax
    
    %Compute the angle quantization depending on the size of the
    %codebook. Note that we quantize the range -1 to 1 uniformly, since
    %this is what impacts the channel directions, and then maps this back
    %onto the user angles.
    if m <= Lmax
        uniformAngles = asin(linspace(-1,1,m));
    else
        uniformAngles = asin(linspace(-1,1,Lmax));
    end

    
    %Prepare to compute different array gains and quantized channel vectors
    arrayGainsTmp = zeros(m,nbrOfRealizations);
    
    %Compute the array gains for all user directions for each of the
    %quantized channel vectors
    for n = 1:length(uniformAngles);
        quantizedChannel = exp( (0:m-1)' * -2i*pi*interAntennaDistance*sin(uniformAngles(n)) ); %Compute the quantized channel direction for a given angle
        arrayGainsTmp(n,:) = abs(quantizedChannel'*H_LoS(1:m,:)).^2 /m; %Compute the corresponding array gain
    end
    
    %Compute the average array gain when each user picks the best channel
    %direction from the codebook
    arrayGainLoS(m) = mean(max(arrayGainsTmp,[],1));
    
end



%%Monte-Carlo simulations for open-loop beamforming with Rayleigh fading

%Generate channel realizations with Rayleigh fading
H_Rayleigh = (randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations))/sqrt(2);

%Prepare to compute different array gains
arrayGainRayleigh = zeros(size(arrayGainMassiveMIMO));


%Go through the range of number of service antennas
for m = 1:Mmax
    
    %Compute the average array gain of the quantized Rayleigh fading, by
    %utilizing the channel distribution is rotationally invariant. This
    %effectively means that quantizing the channel using L orthogonal
    %dimensions reduce the number of antennas to min(L,M).
    if m <= Lmax
        arrayGainRayleigh(m) = mean(max(abs(H_Rayleigh(1:m,:)).^2,[],1),2);
    else
        arrayGainRayleigh(m) = arrayGainRayleigh(m-1);
    end
    
end



%Plot Figure 2b from the paper
figure; hold on; box on;

plot(1:Mmax,arrayGainMassiveMIMO,'r')
plot(1:Mmax,arrayGainLoS,'b-.')

plot(1:Mmax,arrayGainRayleigh,'k--')

xlabel('Number of service antennas (M) ')
ylabel('Average array gain');

plot(Lmax,arrayGainLoS(Lmax),'b*');
plot(Lmax,arrayGainRayleigh(Lmax),'k*');

legend('Massive MIMO: Measured channels','Open-loop: LoS-ULA','Open-loop: Isotropic (Rayleigh)','Location','NorthWest');
