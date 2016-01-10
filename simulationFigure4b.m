%This Matlab script can be used to generate Figure 4b, in the article:
%
%Emil Bjornson, Erik G. Larsson, Thomas L. Marzetta, "Massive MIMO: Ten
%Myths and One Critical Question," IEEE Communications Magazine, To appear.
%
%Download article: http://arxiv.org/pdf/1503.06854
%
%This is version 1.0 (Last edited: 2016-01-10)
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

%Maximal number of Monte-Carlo realizations with random user locations
monteCarloRealizations = 5000;

%Maximal number of users per cell
K = 100;

%Generate range of user numbers
Kvalues = 1:K;

%Select number of BS antennas
M = 100;

%Pathloss exponent
kappa = 3.7;

%SNR value (-5 dB)
SNR = 10^(-5/10);

%Length of coherence block
tau = 200;

%Number of tiers of hexagonals that are simulated, around the desired cell
tiers = 5;

%Percentage of the radius inside the cell where no UEs are allowed
forbiddenRegion = 0.14;

%Define intersite distance in a normalized scale
intersiteDistance = 2;
intersiteDistanceHalf = intersiteDistance/2;

dmax = intersiteDistanceHalf; %Normalized cell radius
dmin = dmax * forbiddenRegion; %Normalized shortest distance from a BS


%%Begin Monte-Carlo simulations

%Placeholders for storing spectral efficiencies from Monte-Carlo
%simulations
averageSEsReuse1_MR = zeros(monteCarloRealizations,length(Kvalues));
averageSEsReuse3_MR = zeros(monteCarloRealizations,length(Kvalues));
averageSEsReuse4_MR = zeros(monteCarloRealizations,length(Kvalues));

averageSEsReuse1_ZF = zeros(monteCarloRealizations,length(Kvalues));
averageSEsReuse3_ZF = zeros(monteCarloRealizations,length(Kvalues));
averageSEsReuse4_ZF = zeros(monteCarloRealizations,length(Kvalues));




%Go through all Monte-Carlo realizations
for n = 1:monteCarloRealizations
    
    %Display simulation progress
    disp(['Realization ' num2str(n) ' out of ' num2str(monteCarloRealizations)]);
    
    
    %Vector with BS positions
    BSpositions = zeros(4,1);
    
    %Vectors that store which set of pilots that each BS uses (when there
    %is non-universal pilot reuse)
    reusePattern3 = zeros(4,1); %Pilot reuse 3
    reusePattern4 = zeros(4,1); %Pilot reuse 4
    
    %BS index of the next BS to be deployed
    itr = 1;
    
    
    %Deploy hexagonal cells in "tiers" number of tiers
    for dim1 = 0:tiers
        for dim2 = 0:tiers-dim1
            
            if (dim1 == 0) || (dim2>=1)
                
                %Compute a BS location
                BSloc = sqrt(3)*dim2*intersiteDistanceHalf*1i + sqrt(3)*dim1*intersiteDistanceHalf*exp(1i*pi*(30/180));
                
                
                %Special: The first BS is placed in the origin (this is
                %where the performance is computed)
                if (dim1 == 0) && (dim2 == 0)
                    BSpositions(itr) = BSloc;
                    reusePattern3(itr) = 1;
                    reusePattern4(itr) = 1;
                    itr = itr+1;
                    
                else
                    
                    %Compute the current BS location
                    basis = [3*intersiteDistanceHalf/2 0; sqrt(3)*intersiteDistanceHalf/2 sqrt(3)*intersiteDistanceHalf];
                    rotation = [cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)];
                    
                    %Deploy BSs in all six directions from the center cell
                    for direction = 1:6
                        
                        %Compute the reuse pattern identity for reuse 3
                        if mod(2*dim1+dim2,3)==0
                            reusePattern3(itr) = 1;
                        elseif mod(mod(2*dim1+dim2,3)+direction,2)==1
                            reusePattern3(itr) = 2;
                        elseif mod(mod(2*dim1+dim2,3)+direction,2)==0
                            reusePattern3(itr) = 3;
                        end
                        
                        %Compute the reuse pattern identity for reuse 4
                        if mod(dim1,2)==0 && mod(dim2,2)==0
                            reusePattern4(itr) = 1;
                        elseif direction == 1
                            reusePattern4(itr) = 1+mod(2*dim1+mod(dim2,2),4);
                        else
                            dims = round(basis\(rotation^(direction-1)*basis*[dim1; dim2]));
                            reusePattern4(itr) = 1+mod(2*dims(1)+mod(dims(2),2),4);
                        end
                        
                        %Deploy a BS
                        BSpositions(itr) = BSloc;
                        itr = itr+1;
                        
                        %Rotate the BS location to consider the next
                        %direction of the six directions from the center
                        %cell
                        BSloc = BSloc .* exp(1i*2*pi*60/360);
                        
                    end
                end
                
            end
            
        end
    end
    
    %Compute the final number of BSs
    nbrBSs = itr-1;
    
    %Prepare to put out UEs in the cells
    UEpositions = zeros(K,nbrBSs);
    
    %Initiate matrices where first and second order interference are computed
    interference1reuse1 = zeros(K,1);
    interference2reuse1 = zeros(K,1);
    interference1reuse3 = zeros(K,3);
    interference2reuse3 = zeros(K,3);
    interference1reuse4 = zeros(K,4);
    interference2reuse4 = zeros(K,4);
    
    
    %Go through all the cells
    for j = 1:nbrBSs
        
        %Generate UE locations randomly with uniform distribution inside the cells
        nbrToGenerate = K; %Number of UE locations left to generate
        notFinished = true(K,1); %Indices of the UE locations that are left to generate
        
        
        %Iterate the generation of UE locations until all of them are inside a
        %hexagonal cell
        while nbrToGenerate>0
            
            %Generate new UE locations uniformly at random in a circle of radius dmax
            UEpositions(notFinished,j) = sqrt( rand(nbrToGenerate,1)*(dmax^2-dmin^2)+ dmin^2 ) .* exp(1i*2*pi*rand(nbrToGenerate,1));
            
            %Check which UEs that are inside a hexagonal and declare as finished
            finished = checkHexagonal(UEpositions(:,j)',dmax);
            
            %Update which UEs that are left to generate
            notFinished = (finished==false);
            
            %Update how many UEs that are left to generate
            nbrToGenerate = sum(notFinished);
            
        end
        
        %Finalize UE locations by translating them around the serving BS
        UEpositions(:,j) = UEpositions(:,j) + BSpositions(j);
        
        %Compute the distance from the users in cell j to BS j
        distancesSquaredBSj = abs(UEpositions(:,j) - BSpositions(j));
        
        
        %Focus on interference caused to the first BS (in the center)
        l = 1;
        
        %Compute the distance from the users in cell j to BS l
        distancesSquaredBSl = abs(UEpositions(:,j) - BSpositions(l));
        
        
        %Compute inteference terms of the types that show up in Eq. (7)
        interference1reuse1(:,l) = interference1reuse1(:,l) + (distancesSquaredBSj./distancesSquaredBSl).^(kappa);
        interference2reuse1(:,l) = interference2reuse1(:,l) + (distancesSquaredBSj./distancesSquaredBSl).^(2*kappa);
        
        interference1reuse3(:,reusePattern3(j)) = interference1reuse3(:,reusePattern3(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(kappa);
        interference2reuse3(:,reusePattern3(j)) = interference2reuse3(:,reusePattern3(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(2*kappa);
        
        
        interference1reuse4(:,reusePattern4(j)) = interference1reuse4(:,reusePattern4(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(kappa);
        interference2reuse4(:,reusePattern4(j)) = interference2reuse4(:,reusePattern4(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(2*kappa);
        
    end
    
    
    %Go through all different number of users
    for kind = 1:length(Kvalues)
        
        k = Kvalues(kind);
        
        %Compute different pilot lengths (depending on pilot reuse factor)
        B = k;
        B3 = 3*k;
        B4 = 4*k;
        
        %Prepare to store the SINRs
        SINRreuse1_MR = zeros(k,1);
        SINRreuse3_MR = zeros(k,1);
        SINRreuse4_MR = zeros(k,1);
        
        SINRreuse1_ZF = zeros(k,1);
        SINRreuse3_ZF = zeros(k,1);
        SINRreuse4_ZF = zeros(k,1);
        
        
        %Compute uplink performance in center cell
        j = 1;
        
        %Compute the SINRs for MR
        SINRreuse1_MR(:,j) = M *ones(k,1) ./ ( (sum(interference1reuse1(1:k,j)) + 1/SNR)*(interference1reuse1(1:k,j)+1/(B*SNR)) + M*(interference2reuse1(1:k,j)-1));
        SINRreuse3_MR(:,j) = M *ones(k,1) ./ ( (sum(interference1reuse1(1:k,j)) + 1/SNR)*(interference1reuse3(1:k,j)+1/(B3*SNR)) + M*(interference2reuse3(1:k,j)-1));
        SINRreuse4_MR(:,j) = M *ones(k,1) ./ ( (sum(interference1reuse1(1:k,j)) + 1/SNR)*(interference1reuse4(1:k,j)+1/(B4*SNR)) + M*(interference2reuse4(1:k,j)-1));
        
        %Compute the SINRs for ZF
        if (M-k)>0
            SINRreuse1_ZF(:,j) = (M-k) *ones(k,1) ./ ( (sum(interference1reuse1(1:k,j)) - sum(interference2reuse1(1:k,j)./(interference1reuse1(1:k,j)+1/(B*SNR))) + 1/SNR)*( interference1reuse1(1:k,j)  +1/(B*SNR)) + (M-k)*(interference2reuse1(1:k,j)-1));
            SINRreuse3_ZF(:,j) = (M-k) *ones(k,1) ./ ( (sum(interference1reuse1(1:k,j)) - sum(interference2reuse3(1:k,j)./(interference1reuse3(1:k,j)+1/(B3*SNR))) + 1/SNR)*( interference1reuse3(1:k,j) +1/(B3*SNR)) + (M-k)*(interference2reuse3(1:k,j)-1));
            SINRreuse4_ZF(:,j) = (M-k) *ones(k,1) ./ ( (sum(interference1reuse1(1:k,j)) - sum(interference2reuse4(1:k,j)./(interference1reuse4(1:k,j)+1/(B4*SNR))) + 1/SNR)*( interference1reuse4(1:k,j) +1/(B4*SNR)) + (M-k)*(interference2reuse4(1:k,j)-1));
        end
        
        
        %Compute the average SEs over the network by Monte-Carlo simulation
        averageSEsReuse1_MR(n,kind) = mean(log2(1+SINRreuse1_MR(:)));
        averageSEsReuse3_MR(n,kind) = mean(log2(1+SINRreuse3_MR(:)));
        averageSEsReuse4_MR(n,kind) = mean(log2(1+SINRreuse4_MR(:)));
        
        averageSEsReuse1_ZF(n,kind) = mean(log2(1+SINRreuse1_ZF(:)));
        averageSEsReuse3_ZF(n,kind) = mean(log2(1+SINRreuse3_ZF(:)));
        averageSEsReuse4_ZF(n,kind) = mean(log2(1+SINRreuse4_ZF(:)));
        
    end
    
    
end


%Plot Figure 4b
figure(4); hold on; box on;


SEreuse1 = Kvalues.*max([(1-Kvalues/tau); zeros(size(Kvalues))]).*mean(averageSEsReuse1_ZF,1);
SEreuse3 = Kvalues.*max([(1-3*Kvalues/tau); zeros(size(Kvalues))]).*mean(averageSEsReuse3_ZF,1);
SEreuse4 = Kvalues.*max([(1-4*Kvalues/tau); zeros(size(Kvalues))]).*mean(averageSEsReuse4_ZF,1);

SEoptimzedZF = max([SEreuse1; SEreuse3; SEreuse4],[],1);

plot(Kvalues,SEoptimzedZF,'k-','LineWidth',1); %Plot ZF


SEreuse1 = Kvalues.*max([(1-Kvalues/tau); zeros(size(Kvalues))]).*mean(averageSEsReuse1_MR,1);
SEreuse3 = Kvalues.*max([(1-3*Kvalues/tau); zeros(size(Kvalues))]).*mean(averageSEsReuse3_MR,1);
SEreuse4 = Kvalues.*max([(1-4*Kvalues/tau); zeros(size(Kvalues))]).*mean(averageSEsReuse4_MR,1);

SEoptimzedMR = max([SEreuse1; SEreuse3; SEreuse4],[],1);

plot(Kvalues,SEoptimzedMR,'b-.','LineWidth',1); %Plot MR


[SEmaxZF,KmaxZF] = max(SEoptimzedZF);
plot(KmaxZF,SEmaxZF,'kd','LineWidth',1);

[SEmaxMR,KmaxMR] = max(SEoptimzedMR);
plot(KmaxMR,SEmaxMR,'bd','LineWidth',1);


axis([0 100 0 45]);
xlabel('Number of users (K)');
ylabel('Spectral efficiency [bit/s/Hz/cell]');

legend('ZF','MR','Location','NorthEast');
