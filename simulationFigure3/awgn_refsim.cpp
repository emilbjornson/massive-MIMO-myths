// Reference simulation with QPSK modulation over an AWGN channel
// using a rate-1/2 LDPC channel code. Symmetry of the LLR
// distribution is exploited and only zeros are transmitted.
//
// This code is based on an IT++ tutorial program, originally written
// by Erik G. Larsson.  The GNU General Public License (GPL) applies.

#include <itpp/itcomm.h>
#include <itpp/itstat.h>

using std::cout;
using std::endl;
using namespace itpp;
using namespace std;

extern int main(int argc, char **argv)
{ 
  if (argc < 2) {
     it_info("Usage: " << argv[0] << " codec_file.it");
     return 1;
  }
  
  cerr << "Shannon limit (unconstrained capacity): 0 dB at rate 1.0 bpcu" << endl;

  // -- simulation parameters --
  const vec rho_dB = "-2:0.25:50";
  const long int Nbitsmax=100*1000*1000;  // maximum number of bits to ever simulate per SNR point
  double FERmin = 0.002;                 // stop simulating a given method if FER<this value
  int Nfers = 100;                       // move to next SNR point after counting this many frame errors

  // -- channel code parameters --  
  LDPC_Generator_Systematic G;
  LDPC_Code C(argv[1], &G);     // read code from file provided on the command line
  C.set_exit_conditions(50);   // 50 decoder iterations 

  const int Nc=C.get_nvar();                 // number of variable nodes (coded bits)
  const int Nu=C.get_nvar()-C.get_ncheck();  // number of uncoded bits

  it_assert(abs(Nc-2*Nu)<25,"this simulation is hardcoded for rate 1/2");
  const int Ns = Nc/2;    // number of symbols (code rate is 1/2)
 
  // -- initialize QPSK channel --
  ND_UQAM chan;
  chan.set_M(Ns,4);  // Ns uses of a scalar channel with QPSK
  ND_UQAM chan1;
  chan1.set_M(1,4);  // one use of a scalar channel with QPSK
 
  // -- run simulation --
  for (int nsnr=0; nsnr<length(rho_dB); nsnr++) {
    cerr << "Processing SNR point " << nsnr << " out of " << length(rho_dB) << endl;
    double rho = pow(10.0,rho_dB(nsnr)/10.0);

    BERC berc;  // counter for coded BER
    BERC bercu; // counter for uncoded BER
    BLERC ferc; // counter for coded FER
    ferc.set_blocksize(Nc);
    
    long int nbits=0;
    long int nframes=0;
    while (nbits<Nbitsmax) {
      nbits += Nu;
      nframes++;

      // generate zero data
      bvec inputbits = zeros_b(Nu);                
      bvec txbits = zeros_b(Nc);    
      it_assert(length(txbits)==Nc,"error");
      it_assert(length(txbits)%2==0,"error");

      cvec yeq = zeros_c(Ns);
      cvec effgain = zeros_c(Ns);

      // symbol loop
      for (int i=0; i<Ns; i++) {
	// generate s(0) -- the symbol of interest
	const cvec s0 = chan1.modulate_bits(txbits.get(2*i,2*i+1));
	const complex<double> s = s0(0);
      
	// generate received data
	const complex<double> w = randn_c();  // The real and
	                                      // imaginary parts of w
	                                      // are independent and
	                                      // have variances equal
	                                      // to 0.5
	const complex<double> y = sqrt(rho)*s + w;
	const complex<double> shat = y; 
	
	yeq(i) = shat;   // equivalent AWGN channel with noise variance 1.0 per complex dimension
	effgain(i) = sqrt(rho);
      }

      // run QPSK demodulator
      QLLRvec llr_apr=zeros_i(Nc);
      QLLRvec llr_d=zeros_i(Nc);
      chan.demodulate_soft_bits(yeq,effgain,1.0,llr_apr,llr_d);
      // The noise variance provided to the ND_UQAM demodulator is the
      // "noise variance per complex dimension, i.e. the sum of real
      // and imaginary parts".

      // decode
      QLLRvec llr_out=zeros_i(Nc);
      int iter_d = C.bp_decode(llr_d,llr_out);
      bvec decoded_bits = llr_out<0;
      bercu.count(txbits,llr_d<0);
      berc.count(txbits,decoded_bits);
      ferc.count(txbits,decoded_bits);

      // check whether it is time to terminate the simulation
      if (ferc.get_errors()>Nfers) { break;}      
    }
    
     cerr << "----------------------" << endl;
     cerr << "rho: " << rho_dB(nsnr) << " dB. Simulated " << nframes << " frames consisting of " << nbits << " bits.";
     cerr << " Uncoded BER: " << bercu.get_errorrate()  
          << " Coded BER: " << berc.get_errorrate() 
          << " FER: " << ferc.get_errorrate() << endl;
     cerr.flush();

     cout << rho_dB(nsnr) << "       " << bercu.get_errorrate()  
	  << "    " << berc.get_errorrate()  
	  << "    "  << ferc.get_errorrate() << endl;
     cout.flush();

     // check whether it is time to stop simulating (stop when reached the min FER of interest)
     if (ferc.get_errorrate()<FERmin)  {  break; } 
  }
  
  return 0;
}
