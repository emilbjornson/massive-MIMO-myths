
F_DEBUG =  -I/home/egl/IT++DEBUG/include -g
F_OPT =  -I/home/egl/IT++DEBUG/include -O3 -DNDEBUG

L_DEBUG = -L/home/egl/IT++DEBUG/lib -litpp_debug /usr/lib/i386-linux-gnu/libfftw3.so /usr/lib/i386-linux-gnu/libm.a  /usr/lib/libblas.so 
L_OPT = -L/home/egl/IT++RELEASE/lib -litpp /usr/lib/i386-linux-gnu/libfftw3.so /usr/lib/i386-linux-gnu/libm.a  /usr/lib/libblas.so

all: mami_zerobits awgn_refsim gen_code

mami_zerobits: mami_zerobits.cpp
	g++ $(F_DEBUG) mami_zerobits.cpp -o mami_zerobits_debug  $(L_DEBUG)
	g++ $(F_OPT) mami_zerobits.cpp -o mami_zerobits_opt  $(L_OPT)

awgn_refsim: awgn_refsim.cpp
	g++ $(F_DEBUG) awgn_refsim.cpp -o awgn_refsim_debug  $(L_DEBUG)
	g++ $(F_OPT) awgn_refsim.cpp -o awgn_refsim_opt  $(L_OPT)

gen_code: gen_code.cpp
	g++ $(F_DEBUG) gen_code.cpp -o gen_code_d  $(L_DEBUG)
	g++ $(F_OPT) gen_code.cpp    -o gen_code_o  $(L_OPT)

mami_sim_all:
	nice -19 ./mami_zerobits_opt RU_1000.it >mami_ru1000 &
	nice -19 ./mami_zerobits_opt RU_10000.it >mami_ru10000 &
	nice -19 ./mami_zerobits_opt RU_100000.it >mami_ru100000 &

awgn_refsim_all:
	nice -19 ./awgn_refsim_opt RU_1000.it >awgn_ru1000 &
	nice -19 ./awgn_refsim_opt RU_10000.it >awgn_ru10000 &
	nice -19 ./awgn_refsim_opt RU_100000.it >awgn_ru100000 &
