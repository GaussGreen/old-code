/* -----------------------------------------------------------
   FILENAME : num_h_spectrunc.h

   PURPOSE	: generate Brownian path through a spectral (PCA)
              truncation
   ----------------------------------------------------------- */	
#ifndef NUM_H_SPECTRUNC_H
#define NUM_H_SPECTRUNC_H


Err SpecTruncCube( 
			double  ***rand,
			double    *time_at_steps,
			double     precision,
			long       first_path, long last_path,
			long       first_brow, long last_brow,
			long       first_step, long last_step,
			long      *seed);

#endif
