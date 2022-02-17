#ifndef Pattern_included

#define Pattern_included


/* A vector in a linked list of vectors */
struct poll_vector {
	double *vector;    /* vector */
	struct poll_vector *next;
};




void pollstep(int n, int pi, double (*objf)(), double *lb, double *ub, struct poll_vector **last_sucess);
void init_D(int);
void print_D(int);
void print_poll_vector(int, double *);
void init_pattern(int);
void clean_pattern();
void clean_D();
double try_poll_point(int, double *, double (*objf)(), double *, double *);

#endif
