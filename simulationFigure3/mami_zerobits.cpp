// Simulation of massive MIMO uplink with QPSK modulation, matched
// filter receiver and MMSE channel estimation and a rate-1/2 LDPC
// channel code. All terminals use orthogonal pilots. The symmetry of
// the LLR distribution is exploited, and only zeros are transmitted
// by the terminal whose data are decoded, circumventing the use of a
// generator matrix.
//
// Copyright Erik G. Larsson, 2014-2016
// Inquiries may be directed to erik.g.larsson@liu.se
// The GNU General Public Licence applies.

#include <itpp/itcomm.h>
#include <itpp/itstat.h>

using std::cout;
using std::endl;
using namespace itpp;
using namespace std;

// -----------
// Find the SNR required to achieve a given rate, as defined by the
// capacity lower bound in Table 3.1, upper right corner, in
// Fundamentals of Massive MIMO by T. Marzetta, E.G. Larsson, H. Yang
// and H.Q. Ngo (Cambridge University Press, 2016)
// -----------
double find_shannon(int M,
		    int K, 
		    double taup,
		    double R_target)
{
  double range_u=30;  // unit is dB
  double range_l=-20; 
  double rho_dB;
  for (int n=0; n<50; n++) {
    rho_dB=(range_l+range_u)/2;
    const double rho=pow(10,rho_dB/10);
    double gamma = taup*rho/(1+taup*rho);
    const double R=log2(1+M*rho*gamma/(1+rho*K));
    if (R>R_target) { range_u=rho_dB; } else {  range_l=rho_dB; }
  }
  
  return rho_dB;
};

// ------------
// Main simulation program
// ------------
extern int main(int argc, char **argv)
{ 
  if (argc < 2) {
     it_info("Usage: " << argv[0] << " codec_file.it");
     return 1;
  }
  
  // -- Massive MIMO parameters --
  const int M=100;       // number of transmitter antennas
  const int K=30;        // number of terminals
  const double taup=K;   // pilot sequence length
  const int Tc=1;        // channel coherence (number of symbols that see the same channel)

  cerr << "Shannon limit: " << find_shannon(M,K,taup,1.0) << " dB at rate 1.0 bpcu" << endl;

  // -- simulation parameters --
  const vec rho_dB = "-16:0.25:50";
  const long int Nbitsmax=100*1000*1000; // maximum number of bits to ever simulate per SNR point
  double FERmin = 0.002;                 // stop simulating a given method if FER<this value
  int Nfers = 100;                       // move to next SNR point after counting this many frame errors

  // -- channel code parameters --  
  LDPC_Generator_Systematic G;
  LDPC_Code C(argv[1], &G);     // read code from file provided on the command line
  C.set_exit_conditions(50);    // 50 decoder iterations 

  // -- initialize --
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

      // generate zero data for the terminal of interest
      bvec inputbits = zeros_b(Nu);                
      bvec txbits = zeros_b(Nc);    
      it_assert(length(txbits)==Nc,"error");
      it_assert(length(txbits)%2==0,"error");

      cvec yeq = zeros_c(Ns);
      cvec effgain = zeros_c(Ns);
      cmat H, Hhat;

      // symbol loop
      for (int i=0; i<Ns; i++) {
	if ((i%Tc)==0) {
	  // generate new channel and channel estimate
	  H = randn_c(M,K); 
	  const cmat Wh = randn_c(M,K);
	  Hhat = sqrt(taup*rho)/(1+taup*rho)*(sqrt(taup*rho)*H + Wh);
	}

	// generate symbols; s(0) is the symbol of interest, all
	// others are interference, constructed from random bits
	cvec s = zeros_c(K);
	s.set_subvector(0,chan1.modulate_bits(txbits.get(2*i,2*i+1)));
	for (int k=1; k<K; k++) {
	  s.set_subvector(k,chan1.modulate_bits(randb(2)));
	}
      
	// generate received data
	const cvec w = randn_c(M);  // The real and imaginary parts of
				    // w are independent and have
				    // variances equal to 0.5
	const cvec y = sqrt(rho)*H*s + w;
	
	// receiver processing (gamma, A, rho, y have the same meaning as in the Massive MIMO book)
	const double gamma = taup*rho / (1+taup*rho);
	const cmat A = Hhat/sqrt(gamma);
	const cvec a0 = A.get_col(0);
	const complex<double> shat = conj(a0) * y / (sqrt(rho*gamma) * sum(sqr(a0)));
	
	// compute effective noise variance, for LLR normalization
	double lasttermdenom=0;
	for (int k=1; k<K; k++) {
	  lasttermdenom += rho*gamma*sqr(conj(a0)*A.get_col(k)) / sum(sqr(a0));
	}
	const double SINR = rho*gamma*sum(sqr(a0)) / (1 + rho*K*(1-gamma) + lasttermdenom);

	yeq(i) = shat * sqrt(SINR);   // equivalent AWGN channel with noise variance 1.0 per complex dimension
	effgain(i)=sqrt(SINR);
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
