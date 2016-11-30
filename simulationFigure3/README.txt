The programs in this directory reproduce the numerical results
reported in Myth 4 in the paper "Massive MIMO: Ten Myths and One
Critical Question" by E. Bjornson, E.G. Larsson and T. Marzetta (IEEE
Communications Magazine, 2015).

Inquiries may be directed to Erik G. Larsson, erik.g.larsson@liu.se

The instructions herein assume a Ubuntu 16.04.1 LTS Linux
environment. Other platforms or Linux distributions may work too, but
tweaking of the linking flags in the Makefile may be necessary.

The steps required to reproduce the results are as follows:

1. Make sure the following packages are installed:

liblapack3 liblapack-dev libfftw3-dev libfftw3-bin cmake

2. Install IT++ (itpp.sf.net), release 4.3.1.  (The code may work with
other releases, but this has not been tested.) The steps are as follows.

- Download googletest-master.zip into ~, unzip, cmake . and make
- Download it++ 4.3.1 into ~, and unpack
- cd itpp-4.3.1
- mkdir build
- cd build
- rm -rf *
- cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/egl/IT++RELEASE   -DGTEST_DIR=/home/egl/googletest-master/googletest
- nice -19 make && make install
- rm -rf *
- cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/home/egl/IT++DEBUG   -DGTEST_DIR=/home/egl/googletest-master/googletest
- nice -19 make && make install
- cd gtests
- ./itpp_gtests
- Verify that all tests passed

3. If IT++ is installed as described above, the Makefile should work
as is. Otherwise, adjust the linking paths therein as necessary.

4. Run make all

5. Set the linker path,

export LD_LIBRARY_PATH=/home/egl/IT++DEBUG/lib:/home/egl/IT++RELEASE/lib/

6. Run "gen_code_opt" to generate the LDPC codes. This takes a
while. The computer must have at least 2 GB of memory.

7. Run "make awgn_refsim_all" to run the AWGN reference
simulations. This may take overnight. The purpose is to verify the
correct operation of the codec. Run "gnuplot berplot_awgn.gp" to
generate the file "awgn.eps", with AWGN performance of the code. Check
that this figure looks reasonable. 

8. Run "make mami_zerobits_opt" to run the Massive MIMO
simulation. This may take overnight, even on a top-of-the-line
PC. Then run "gnuplot berplot_mami.gp" to generate the figure in the
paper. The C++ program also computes the SNR required to meet a given
value of the capacity lower bound.

Further comments can be found inside of the .cpp files.
