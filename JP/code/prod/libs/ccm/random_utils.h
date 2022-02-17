/* $Header$ */
int CreateGaussianRandomSequence(
    double *Sequence,           // (O) Ptr to object pointer        
	long    aSeed,              // (I) Seed for gasdev and perm of partitions   
    int     numberDimensions,   // (I) Number of dimensions                     
	int     numberPaths);       // (I) Number of paths to extract               

int  CreateSobolSequence(
    double  *Sequence,          // (O) Ptr to object pointer                   
	long    aSeed,              // (I) Seed for gasdev and perm of partitions   
    int     numberDimensions,   // (I) Number of dimensions                     
	int     numberPaths);       // (I) Number of paths to extract               


int  CreateUniformSequence(
    double  *Sequence,          // (O) Ptr to object pointer                   
	long    aSeed,              // (I) Seed for gasdev and perm of partitions   
    int     numberDimensions,   // (I) Number of dimensions                     
	int     numberPaths);       // (I) Number of paths to extract               

void *CreateRandomGenerator(long aSeed);
void RandomGeneratorFree(void *A);
double RandomGeneratorGet(void *A);
