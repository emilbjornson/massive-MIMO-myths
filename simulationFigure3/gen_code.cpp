// Generate some LDPC codes optimized for the AWGN channel.  Generator
// matrices are not computed.

// This code is adapted from an IT++ tutorial program, originally
// written by Erik G. Larsson.  The GNU General Public License (GPL)
// applies.

#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;

int main(int argc, char **argv)
{
  { // This generates a random regular (3,6) code with 2000 bits
    LDPC_Parity_Regular H;
    H.generate(2000, 3, 6,
	       "rand",  // random unstructured matrix
	       "2000 10");   // optimize girth
    H.cycle_removal_MGW(8);
    H.display_stats();
    LDPC_Generator_Systematic G1(&H);
    LDPC_Code C1(&H,&G1);
    C1.save_code("36-reg-code-2000bits.it");
  }

  // Irregular 1/2-rate codes optimized for the AWGN channel. The
  // degree distributions are taken from Richardson & Urbanke,
  // Trans. IT 2001.

  { // 1000 bits
    cout << "========= IRREGULAR CODE 1000 BITS ==========" << endl;
    LDPC_Parity_Irregular H;
    H.generate(1000,
               "0 0.27684 0.28342 0 0 0 0 0 0.43974",
               "0 0 0 0 0 0.01568 0.85244 0.13188",
               "rand",   // random unstructured matrix
               "500 8"); // optimize girth
    LDPC_Code C(&H);
    C.save_code("RU_1000.it");
  }

  {  // 10000 bits (takes a few minutes to run)
    cout << "========= IRREGULAR CODE 10000 BITS ==========" << endl;
    LDPC_Parity_Irregular H;
    H.generate(10000,
               "0 0.21991 0.23328 0.02058 0 0.08543 0.06540 0.04767 0.01912 "
               "0 0 0 0 0 0 0 0 0 0.08064 0.22798",
               "0 0 0 0 0 0 0 0.64854 0.34747 0.00399",
               "rand",  // random unstructured matrix
               "150 8"); // optimize
    LDPC_Code C(&H);
    C.save_code("RU_10000.it");
  }

  { // 100000 bits (takes a while to run)
    cout << "========= IRREGULAR CODE 100000 BITS ==========" << endl;
    LDPC_Parity_Irregular H;
    H.generate(100000,
               "0 0.1712 0.21053 0.00273 0 0 0.00009 0.15269 0.09227 "
               "0.02802 0 0 0 0 0.01206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
               "0.07212 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.25830",
               "0 0 0 0 0 0 0 0 0.33620 0.08883 0.57497",
               "rand",
               "40 4");  // less aggressive optimization
    LDPC_Code C(&H);
    C.save_code("RU_100000.it");
  }


  { // 1000000 bits (THIS CODE REQUIRES ABOUT 450 MB TO STORE AND 2GB
    // INTERNAL MEMORY TO GENERATE)
    cout << "========= IRREGULAR CODE 1000000 BITS ==========" << endl;
    LDPC_Parity_Irregular H;
    H.generate(1000000,
               "0 0.1712 0.21053 0.00273 0 0 0.00009 0.15269 0.09227 "
               "0.02802 0 0 0 0 0.01206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
               "0.07212 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.25830",
               "0 0 0 0 0 0 0 0 0.33620 0.08883 0.57497",
               "rand",
               "0 0");  // no optimization here
    LDPC_Code C(&H);
    C.save_code("RU_1000000.it");
  }

  return 0;

}
