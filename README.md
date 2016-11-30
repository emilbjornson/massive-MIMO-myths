Massive MIMO: Ten Myths and One Critical Question
==================

This is a code package is related to the follow scientific magazine article:

Emil Björnson, Erik G. Larsson, Thomas L. Marzetta, “[Massive MIMO: Ten Myths and One Critical Question](http://arxiv.org/pdf/1503.06854),” IEEE Communications Magazine, vol. 54, no. 2, pp. 114-123, February 2016.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


##Abstract of Article

Wireless communications is one of the most successful technologies in modern years, given that an exponential growth rate in wireless traffic has been sustained for over a century (known as Cooper’s law). This trend will certainly continue driven by new innovative applications; for example, augmented reality and internet-of-things.

Massive MIMO (multiple-input multiple-output) has been identified as a key technology to handle orders of magnitude more data traffic. Despite the attention it is receiving from the communication community, we have personally witnessed that Massive MIMO is subject to several widespread misunderstandings, as epitomized by following (fictional) abstract:

“<i>The Massive MIMO technology uses a nearly infinite number of high-quality antennas at the base stations. By having at least an order of magnitude more antennas than active terminals, one can exploit asymptotic behaviors that some special kinds of wireless channels have. This technology looks great at first sight, but unfortunately the signal processing complexity is off the charts and the antenna arrays would be so huge that it can only be implemented in millimeter wave bands.</i>”

The statements above are, in fact, completely false. In this overview article, we identify ten myths and explain why they are not true. We also ask a question that is critical for the practical adoption of the technology and which will require intense future research activities to answer properly. We provide references to key technical papers that support our claims, while a further list of related overview and technical papers can be found at the Massive MIMO Info Point: [http://massivemimo.eu](http://massivemimo.eu)


##Content of Code Package

This package contains C++ code that generates Figure 3 (Myth 4), see the folder simulationFigure3 for details.
This package also contains Matlab-scripts that generate 5 of the simulation figures in the article: simulationFigure2a.m, simulationFigure2b.m, simulationFigure4a.m, simulationFigure4b.m, and simulationFigure5.m

The package contains 3 additional Matlab functions: checkHexagonal.m, function_capacity_broadcast.m and functionHeuristicPowerAllocation.m. These functions are called by the Matlab scripts.

The convex optimization problems are implemented using the modeling language [CVX](http://cvxr.com/cvx/).

See each file for further documentation.


##Acknowledgements

The writing of this article was supported by the EU FP7 under ICT-619086 (MAMMOET) and by ELLIIT and CENIIT.


##License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
