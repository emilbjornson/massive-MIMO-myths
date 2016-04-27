%This Matlab script can be used to generate Figure 2a, in the article:
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
M = 100; %Number of service antennas
K = 12; %Number of terminals


%%Generate channel realizations

%Generate channel realizations with Rayleigh fading
H_Rayleigh = sqrt(1/2)*(randn(M,K,nbrOfRealizations)+1i*randn(M,K,nbrOfRealizations));


%Generate uniformly distributed random user terminal angles for LoS
%propagation; note that -pi/2 to pi/2 is analytically equivalent to -pi to
%pi since a ULA cannot separate all angles
userAngles = pi*rand(1,K,nbrOfRealizations)-pi/2;

%Define distance between antennas in the ULA, measured in wavelengths
interAntennaDistance = 0.5;

%Generate channel realizations with LoS propagation
H_LoS =  exp( repmat((0:M-1)',[1 K nbrOfRealizations]) .* repmat(-2i*pi*interAntennaDistance*sin(userAngles),[M 1 1]) );


%Placeholders for storing simulation results
capacity_Rayleigh = zeros(nbrOfRealizations,1); %Capacity with K terminals and Rayleigh fading
capacity_Rayleigh_drop = zeros(nbrOfRealizations,1); %Capacity with K-2 terminals and Rayleigh fading

capacity_LoS = zeros(nbrOfRealizations,1); %Capacity with K terminals and LoS propagation
capacity_LoS_drop = zeros(nbrOfRealizations,1); %Capacity with K-2 terminals and LoS propagation


%%Compute simulation results


%Go through all Monte-Carlo realizations
for n = 1:nbrOfRealizations
    
    %Output the progress of the simulation
    if mod(n,10000)==0
        disp(['Realization ' num2str(n) ' out of ' num2str(nbrOfRealizations)]);
    end
    
     %Normalize the Rayleigh fading realizations to achieve isotropically
     %distributed channel directions. The reason for this is to make the
     %favorable propagation curve representative for both channel types.
     H_Rayleigh_norm = sqrt(M)*H_Rayleigh(:,:,n)./repmat(sqrt(sum(abs(H_Rayleigh(:,:,n)).^2,1)),[M 1]);
    
     %Compute the uplink sum capacity for both channel types, using (10.15)
     %in "Fundamentals of Wireless Communication" by David Tse and Pramod
     %Viswanath, under the assumption that both terminals operate at full
     %power to achieve the prescribed SNR.
     capacity_Rayleigh(n) = abs(log2(det(eye(K) + SNR*(H_Rayleigh_norm'*H_Rayleigh_norm))));
     capacity_LoS(n) = abs(log2(det(eye(K) + SNR*(H_LoS(:,:,n)'*H_LoS(:,:,n)))));
     
     %Go through all possible ways of dropping two terminals 
     for k1 = 1:K-1
         
         for k2 = k1+1:K
             
             %Indicies of the active terminals
             index = [1:k1-1 k1+1:k2-1 k2+1:K];
             
             %Compute the capacity in the same way as above, but with only
             %K-2 terminals
             rate_Rayleigh_drop = abs(log2(det(eye(K-2) + SNR*(H_Rayleigh_norm(:,index)'*H_Rayleigh_norm(:,index)))));
             rate_LoS_drop = abs(log2(det(eye(K-2) + SNR*(H_LoS(:,index,n)'*H_LoS(:,index,n)))));
             
             %Store the best result with Rayleigh fading
             if rate_Rayleigh_drop>capacity_Rayleigh_drop(n)
                 capacity_Rayleigh_drop(n) = rate_Rayleigh_drop;
             end
             
             %Store the best result with LoS propagation
             if rate_LoS_drop>capacity_LoS_drop(n)
                 capacity_LoS_drop(n) = rate_LoS_drop;
             end
             
         end
         
     end
     
end

%Generate a vertical line to be used for the favorable propagation case
verticalLine = linspace(0,1,1000);


%Plot Figure 2a from the paper
figure; hold on; box on;

plot(sort(capacity_LoS),linspace(0,1,nbrOfRealizations),'b-.','LineWidth',1);
plot(sort(capacity_Rayleigh),linspace(0,1,nbrOfRealizations),'r','LineWidth',1);
plot(K*log2(1+SNR*M)*ones(size(verticalLine)),verticalLine,'k--','LineWidth',1);

plot(sort(capacity_LoS_drop),linspace(0,1,nbrOfRealizations),'b-.','LineWidth',1);
plot(sort(capacity_Rayleigh_drop),linspace(0,1,nbrOfRealizations),'r','LineWidth',1);
plot((K-2)*log2(1+SNR*M)*ones(size(verticalLine)),verticalLine,'k--','LineWidth',1);

set(gca,'Yscale','log');
xlabel('Sum capacity [bit/symbol]');
ylabel('Cumulative probability');
ylim([0.01 1]);
xlim([42 62]);

legend('LoS-ULA','Isotropic (Rayleigh)','FP channel','Location','NorthWest');
