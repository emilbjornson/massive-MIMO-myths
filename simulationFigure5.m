%This Matlab script can be used to generate Figure 5, in the article:
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

%Define system dimensions
M = logspace(1,3,100); %Range of number of service antennas
antennasPerTerminal = 5; %Fix the ratio M/K
K = M/antennasPerTerminal; %Compute range of number of terminals
tauValues = logspace(2, 4, 100); %Range of coherence block length

%Define OFDM spectral parameters
N = 1200; %Number of subcarriers
subcarrierBandwidth = 15e3; %Bandwidth per subcarrier [Hz]
Btotal = N*subcarrierBandwidth; %Compute effective bandwidth [Hz]

%Define computational efficiency
L_BS = 12.8e9; %flop/W


%%Compute computational complexity

%%Compute FFT complexity (flops)
FFTsize = 1.7*N; %Size of FFT with an oversampling factor

%Number of operations per second with a 4 log2(S) S complexity algorithm,
%where S is the size of the FFT
fftOperations = 4*FFTsize*log2(FFTsize) * subcarrierBandwidth; 


%Placeholders for storing simulation results
complexityZF = zeros(length(tauValues),length(M));
complexityMR = zeros(length(tauValues),length(M));


%Go through different length of the coherence block
for ind = 1:length(tauValues)

    %Extract current coherence block length
    tau = tauValues(ind);

    %Compute the power that scales linearly with M
    perM = fftOperations/L_BS;

    
    %Compute the computational complexity of baseband processing based on
    %the model in "Optimal Design of Energy-Efficient Multi-User MIMO
    %Systems: Is Massive MIMO the Answer?" by Emil Bj?rnson, Luca
    %Sanguinetti, Jakob Hoydis, M?rouane Debbah
    
    %Compute the power that scales as K^3
    perK3_ZF = Btotal/(3*tau*L_BS); %Channel inversion using Cholesky factorization costs K^3/3
    perK3_MR = 0; %No channel inversion in MR
    
    %Compute the power that scales as M K^2
    perMK2_ZF = Btotal*(3)/(tau*L_BS); %Channel inversion using Cholesky factorization costs 3 M K^2
    perMK2_MR = 0; %No channel inversion in MR
    
    %Compute the power that scales as M K
    perMK_ZF = Btotal*(2+1/tau)/L_BS; %Complexity of data transmission and channel estimation
    perMK_MR = Btotal*(2+3/tau)/L_BS; %Complexity of MR precoding computation, data transmission and channel estimation
    
    
    %Compute the power consumption with ZF and MR
    complexityZF(ind,:) = perK3_ZF*K.^3 + perM*M + perMK_ZF*M.*K + perMK2_ZF*M.*K.^2;
    complexityMR(ind,:) = perK3_MR*K.^3 + perM*M + perMK_MR*M.*K + perMK2_MR*M.*K.^2;

    %Remove the points that are not feasible due to impractical overhead
    complexityZF(ind,K>tau) = NaN;
    complexityMR(ind,K>tau) = NaN;
    
end


%Select number of colors in the contour plots
nbrOfContours = 25;


%Plot the MR part of Figure 5
figure(1);

contourf(M,tauValues,10*log10(complexityMR),nbrOfContours)
set(gca,'YScale','log');
set(gca,'XScale','log');

xlabel('Number of service antennas (M)');
ylabel('Coherence block length \tau (symbols)');
title('MR');

%Add a color bar that shows the power consumption in dbW
hcb = colorbar;
colormap(flipud(hot));
caxis([0,30]);
set(gca,'XTickLabel',{'10','100','1000'});
set(gca,'YTickLabel',{'100','1000','10000'});


%Plot the ZF/MMSE part of Figure 5
figure(2);

contourf(M,tauValues,10*log10(complexityZF),nbrOfContours)
set(gca,'YScale','log');
set(gca,'XScale','log');

xlabel('Number of service antennas (M)');
ylabel('Coherence block length \tau (symbols)');
title('ZF/MMSE');

%Add a color bar that shows the power consumption in dbW
colorbar('WestOutside');
colormap(flipud(hot));
caxis([0,30]);
set(gca,'XTickLabel',{'10','100','1000'});
set(gca,'YTickLabel',{'100','1000','10000'});
