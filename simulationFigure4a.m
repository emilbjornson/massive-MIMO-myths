%This Matlab script can be used to generate Figure 4a, in the article:
%
%Emil Bj?rnson, Erik G. Larsson, Thomas L. Marzetta, "Massive MIMO: Ten
%Myths and One Critical Question," IEEE Communications Magazine, To appear.
%
%Download article: http://arxiv.org/pdf/1503.06854
%
%This is version 1.0 (Last edited: 2015-08-29)
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

%Should the sum capacity be calculated? (true /false)
%(Note that it is the calculation of the sum capacity is the main source of
%computational complexity, since we need to numerically optimize the power
%allocation. It can take hours or days to get a sufficient number of
%realizations, while it only takes a few seconds or minutes for heuristic
%precoding and power allocation. This shows clearly why heuristic schemes
%such as ZF and MR are useful in practice.
computeCapacity = true; %true;


%Select the number of realizations in the Monte-Carlo simulations
nbrOfMonteCarloRealizations = 1;

%Define system dimensions
Mrange = 20:10:100; %Range of number of service antennas
K = 20; %Number of terminals

%Define SNR value
PdB = -5; %Power per user, in dB scale
Psum = K*10.^(PdB/10); %Total transmit power, in linear scale




%User weights for (unweighted) sum rate computation
weights = ones(K,1);


%Placeholders for storing simulation results
sumRateMR = zeros(nbrOfMonteCarloRealizations,length(Mrange));
sumRateZF = zeros(nbrOfMonteCarloRealizations,length(Mrange));
sumRateFP = zeros(nbrOfMonteCarloRealizations,length(Mrange));

if computeCapacity == true
    sumCapacity = zeros(nbrOfMonteCarloRealizations,length(Mrange));
end



%Go through the different number of service antennas
for m = 1:length(Mrange)
    
    %Extract the current number of service antennas
    M = Mrange(m);
    
    %Generate of Rayleigh fading channel realizations
    Hall = (randn(K,M,nbrOfMonteCarloRealizations)+1i*randn(K,M,nbrOfMonteCarloRealizations))/sqrt(2);
    
    %Output the simulation progress
    disp(['Progress: M = ' num2str(M)]);
    
    %Go through all channel realizations
    for n = 1:nbrOfMonteCarloRealizations
        
        %Extract the current channel realization
        H = Hall(:,:,n);
        
        %Compute normalized MR precoding
        wMR = H'./repmat(sqrt(sum(abs(H').^2,1)),[M 1]);
        
        %Compute normalized ZF precoding
        wZF = H'/(H*H');
        wZF = wZF./repmat(sqrt(sum(abs(wZF).^2,1)),[M 1]);
        
        
        %Compute a sum rate with MR precoding

        %Optimize the power allocation with MR precoding using Theorem 3.5
        %in "Optimal Resource Allocation in Coordinated Multi-Cell Systems"
        %by Emil Bj?rnson and Eduard Jorswieck
        rhos = diag(abs(H*wMR).^2)';
        powerAllocationMR = functionHeuristicPowerAllocation(rhos,Psum,ones(K,1));
        
        %Calculate sum rate with MR precoding
        W_MR = kron(sqrt(powerAllocationMR),ones(M,1)).*wMR;
        channelGains = abs(H*W_MR).^2;
        signalGains = diag(channelGains);
        interferenceGains = sum(channelGains,2)-signalGains;
        rates = log2(1+signalGains./(interferenceGains+1));
        sumRateMR(n,m) = sum(rates);
        
        
        %Compute a sum rate/capacity with favorable propagation, where MR
        %precoding is optimal and the interference becomes zero
        sumRateFP(n,m) = sum(log2(1+signalGains));
        
        
        %Optimize the power allocation with ZF precoding using Theorem 3.5
        %in "Optimal Resource Allocation in Coordinated Multi-Cell Systems"
        %by Emil Bj?rnson and Eduard Jorswieck
        rhos = diag(abs(H*wZF).^2)';
        powerAllocationwZF = functionHeuristicPowerAllocation(rhos,Psum,ones(K,1));
        
        %Calculate sum rate with ZF precoding
        W_ZF = kron(sqrt(powerAllocationwZF),ones(M,1)).*wZF;
        channelGains = abs(H*W_ZF).^2;
        signalGains = diag(channelGains);
        interferenceGains = sum(channelGains,2)-signalGains;
        rates = log2(1+signalGains./(interferenceGains+1));
        sumRateZF(n,m) = sum(rates);
        

        %Check if the sum capacity should be computed
        if computeCapacity == true
            
            %Compute the sum capacity based on Theorem in "Sum Capacity of
            %the Vector Gaussian Broadcast Channel and Uplink?Downlink
            %Duality" by Pramod Viswanath and David Tse
            sumCapacity(n,m) = function_capacity_broadcast(H,Psum);
            
        end
        
    end
    
end


%Plot Figure 4a from the paper
figure; hold on; box on;

plot(Mrange,mean(sumRateFP,1),'k--','LineWidth',1);

if computeCapacity==true
    plot(Mrange,mean(sumCapacity,1),'ro-','LineWidth',1);
end

plot(Mrange,mean(sumRateZF,1),'b*--','LineWidth',1);
plot(Mrange,mean(sumRateMR,1),'kv-.','LineWidth',1);

if computeCapacity == true
    legend('FP (no interference)','Sum capacity','ZF','MR','Location','NorthWest');
elseif computeCapacity == false
    legend('FP (no interference)','ZF','MR','Location','NorthWest');
end

xlabel('Number of service antennas (M)');
ylabel('Spectral efficiency [bit/s/Hz/cell]');
