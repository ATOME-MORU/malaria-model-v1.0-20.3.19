

#define PI 3.1415927
#define MUTMAX 200
#define NPAR_PER_DEME_FIX 150
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long g_seed;

float inline ran2()
{
	int j; long k;
	static long g_seed2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (g_seed <= 0) {
		if (-(g_seed) < 1) g_seed=1;
		else g_seed = -(g_seed);
		g_seed2=(g_seed);
		for (j=NTAB+7; j>=0; j--) {
			k=(g_seed)/IQ1;
			g_seed=IA1*(g_seed-k*IQ1)-k*IR1;
			if (g_seed < 0) g_seed += IM1;
			if (j < NTAB) iv[j] = g_seed;

		}
		iy=iv[0];
	}
	k=(g_seed)/IQ1;
	g_seed=IA1*(g_seed-k*IQ1)-k*IR1;
	if (g_seed < 0) g_seed += IM1;
	k=g_seed2/IQ2;
	g_seed2=IA2*(g_seed2-k*IQ2)-k*IR2;
	if (g_seed2 < 0) g_seed2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-g_seed2;
	iv[j]=g_seed;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

