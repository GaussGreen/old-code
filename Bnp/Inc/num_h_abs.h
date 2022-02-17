/* -------------------------------------------------------------------------
   FILENAME		: num_h_abs.h

   PURPOSE		: generate randon numbers that are uniformely
   distributed , but with an antithetic noise
   ------------------------------------------------------------------------- */
#ifndef NUM_H_ABS_H
#define NUM_H_ABS_H

Err ABSCube(double ***rand, long first_path, long last_path,
            long first_brownian, long last_brownian, long first_step,
            long last_step, long *seed);

#endif
