#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CSVLib.h>
#include <math.h>
#include <assert.h>

#define MAXLOG 7.09782712893383996732E2
#define Z_MAX 6.0
#define Z_EPSILON 0.000001
#define MIN_CHUNK 64

#ifdef _WIN32
	#define KNRM  " "
	#define KRED  " "
	#define KGRN  " "
	#define KYEL  " "
	#define KBLU  " "
	#define KMAG  " "
	#define KCYN  " "
	#define KWHT  " "
#else
	#define KNRM  "\x1B[0m"
	#define KRED  "\x1B[31m"
	#define KGRN  "\x1B[32m"
	#define KYEL  "\x1B[33m"
	#define KBLU  "\x1B[34m"
	#define KMAG  "\x1B[35m"
	#define KCYN  "\x1B[36m"
	#define KWHT  "\x1B[37m"
#endif

struct pair
{
	double min;
	double max;
};

struct cmdparams
{
	int paired;
	int filter;
	double clevel;
	double pmean;
	char *fname1;
	char *fname2;
};

void array_sort(double *array , int n)
{
	int i,j,temp;

	for(i=0 ; i<n ; i++)
	{
		for(j=0 ; j<n-1 ; j++)
		{
			if(array[j]>array[j+1])
			{
				temp = array[j];
				array[j] = array[j+1];
				array[j+1] = temp;
			}
		}
	}
}

double PerDiff (double f, double s)
{
	if (f == s)
		return 0.0;

	if (f < s)
		return (s - f) / f * 100.0;
	else
		return (f - s) / f * 100.0;
}

double avg (double *buffer, int count)
{
	int i;
	double total = 0;

	for (i = 0; i < count; i++)
		total += buffer[i];

	return total / count;
}

double SDPop (double *buffer, int count)
{
	double mean, standardDeviation = 0.0;
	int i;
	mean = avg (buffer, count);

	for (i = 0; i < count; ++i)
		standardDeviation += pow (buffer[i] - mean, 2);

	return sqrt (standardDeviation / count);
}

double SDSamp (double *buffer, int count)
{
	double mean, standardDeviation = 0.0;
	int i;
	mean = avg (buffer, count);

	for (i = 0; i < count; ++i)
		standardDeviation += pow (buffer[i] - mean, 2);

	return sqrt (standardDeviation / (count - 1));
}

double StandardError(double* buffer, int count)
{
	return SDSamp(buffer, count) / sqrt(count);
}

int SDofDifferences(double *rtn, double* array1, double* array2, int count)
{
	int i;
	double *diff;

	diff = calloc(count, sizeof(double));
	if (!diff) return -1;

	for (i=0; i<count; i++)
	{
		diff[i] = array2[i] - array1[i];
	}

	*rtn = SDSamp(diff,count);
	free (diff);

	return 0;
}

int SEofDifferences(double *rtn,double *array1, double *array2, int count)
{
	int i;
	double *diff;

	diff = calloc(count, sizeof(double));
	if (!diff) return -1;

	for (i=0; i<count; i++)
	{
		diff[i] = array2[i] - array1[i];
	}

	*rtn = SDSamp(diff, count) / sqrt(count);
	free (diff);

	return 0;
}

int MeanofDifferences(double *rtn, double *array1, double *array2, int count)
{
	int i;
	double *diff;

	diff = calloc(count, sizeof(double));
	if (!diff) return -1;

	for (i = 0; i < count; i++)
	{
		diff[i] = array2[i] - array1[i];
	}

	*rtn = avg(diff, count);

	free(diff);

	return 0;
}

double median(double array[], int n)
{
    double median=0;
    
    if(n%2 == 0)
        median = (array[(n-1)/2] + array[n/2])/2.0;
    else
        median = array[n/2];
    
    return median;
}

double R(double *x, double *y, int array_size)
{
	double Mx, My;
	double XMxSum = 0, YMySum = 0;
	double XMxYMySum = 0;
	int i;

	Mx = avg(x, array_size);
	My = avg(y, array_size);

	for (i=0;i<array_size;i++)
	{
		XMxSum += (x[i] - Mx) * (x[i] - Mx);
		YMySum += (y[i] - My) * (y[i] - My);
		XMxYMySum += (x[i] - Mx) * (y[i] - My);
	}

	return XMxYMySum / sqrt(XMxSum * YMySum);
}

double poz (double z)
{
	double  y, x, w;

	if (z == 0.0)
		x = 0.0;
	else
	{
		y = 0.5 * fabs (z);

		if (y >= (Z_MAX * 0.5))
			x = 1.0;
		else if (y < 1.0)
		{
			w = y * y;
			x = ((((((((0.000124818987 * w
			            - 0.001075204047) * w + 0.005198775019) * w
			          - 0.019198292004) * w + 0.059054035642) * w
			        - 0.151968751364) * w + 0.319152932694) * w
			      - 0.531923007300) * w + 0.797884560593) * y * 2.0;
		}
		else
		{
			y -= 2.0;
			x = (((((((((((((-0.000045255659 * y
			                 + 0.000152529290) * y - 0.000019538132) * y
			               - 0.000676904986) * y + 0.001390604284) * y
			             - 0.000794620820) * y - 0.002034254874) * y
			           + 0.006549791214) * y - 0.010557625006) * y
			         + 0.011630447319) * y - 0.009279453341) * y
			       + 0.005353579108) * y - 0.002141268741) * y
			     + 0.000535310849) * y + 0.999936657524;
		}
	}

	return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}

double critz (double p)
{
	double  minz = -Z_MAX;
	double  maxz = Z_MAX;
	double  zval = 0.0;
	double  poz (), pval;

	if (p <= 0.0 || p >= 1.0)
		return (0.0);

	while (maxz - minz > Z_EPSILON)
	{
		pval = poz (zval);

		if (pval > p)
			maxz = zval;
		else
			minz = zval;

		zval = (maxz + minz) * 0.5;
	}

	return (fabs (zval));
}

int RemoveOutliersUnpaired (double *in, double **out, int n, double sensitivity)
{
	int c = 0;
	int i;
	double test;
	double a = avg (in, n);
	double std = SDSamp (in, n);
	*out = malloc (n * sizeof (double));

	if (*out == NULL)
		return 0;

	for (i = 0; i < n; i++)
	{
		test = n * erfc (fabs (in[i] - a) / std);

		if (test < sensitivity)
			printf ("\n%sThrew out #%d -> %s%.6g", KGRN, i + 1, KYEL, in[i]);
		else
			(*out)[c++] = in[i];
	}

	return c;
}

int RemoveOutliersPaired(double *inp1, double **out1, double *inp2, double **out2, int n, double sensitivity)
{
		int c = 0;
		int i;
		double testa,testb;
		double a = avg(inp1, n);
		double stda = SDSamp(inp1, n);
		double b = avg(inp2, n);
		double stdb = SDSamp(inp2, n);
		*out1 = malloc (n * sizeof (double));
		*out2 = malloc (n * sizeof (double));

		for (i = 0; i < n; i++)
		{
				testa = n * erfc(fabs(inp1[i] - a) / stda);
				testb = n * erfc(fabs(inp2[i] - b) / stdb);

				if (testa < sensitivity || testb < sensitivity)
				{
					printf ("\n%sThrew out pair #%d -> %s%.6g , %.6g", KGRN, i + 1, KYEL, inp1[i], inp2[i]);
				}
				else
				{
					(*out1)[c] = inp1[i];
					(*out2)[c] = inp2[i];
					c++;
				}
		}

		return c;
}

double PfromT (const double welch_t_statistic, const double dof)
{
	const double a = dof / 2;
	double value = dof / (welch_t_statistic * welch_t_statistic + dof);

	if ((isinf (value) != 0) || (isnan (value) != 0)) return 1.0;

	printf ("\n%sT = %s%.6g", KGRN, KYEL, welch_t_statistic);
	printf ("\n%sdf = %s%.6g\n", KGRN, KYEL, dof);

	const double beta = lgammal (a) + 0.57236494292470009 - lgammal (a + 0.5);
	const double acu = 0.1E-14;
	double ai;
	double cx;
	int indx;
	int ns;
	double pp;
	double psq;
	double qq;
	double rx;
	double temp;
	double term;
	double xx;

	if (value < 0.0 || 1.0 < value)
		return value;

	if (value == 0.0 || value == 1.0)
		return value;

	psq = a + 0.5;
	cx = 1.0 - value;

	if (a < psq * value)
	{
		xx = cx;
		cx = value;
		pp = 0.5;
		qq = a;
		indx = 1;
	}
	else
	{
		xx = value;
		pp = a;
		qq = 0.5;
		indx = 0;
	}

	term = 1.0;
	ai = 1.0;
	value = 1.0;
	ns = (int) (qq + cx * psq);
	rx = xx / cx;
	temp = qq - ai;

	if (ns == 0)
		rx = xx;

	for (;;)
	{
		term = term * temp * rx / (pp + ai);
		value = value + term;;
		temp = fabs (term);

		if (temp <= acu && temp <= acu * value)
		{
			value = value * exp (pp * log (xx) + (qq - 1.0) * log (cx) - beta) / pp;

			if (indx)
				value = 1.0 - value;

			break;
		}

		ai = ai + 1.0;
		ns = ns - 1;

		if (0 <= ns)
		{
			temp = qq - ai;

			if (ns == 0)
				rx = xx;
		}
		else
		{
			temp = psq;
			psq = psq + 1.0;
		}
	}

	return value;
}

double PValue(double *array1, int array_size, double u)
{
	double mean;
	double std;
	double welch_t_statistic;
	double dof = array_size - 1;

	mean = avg(array1, array_size);
	std = SDPop(array1, array_size);
	welch_t_statistic = (mean - u) / (std / sqrt(array_size - 1));

	return PfromT(welch_t_statistic, dof);
}

double PValueUnpaired (double *array1, const size_t array1_size, double *array2, const size_t array2_size)
{
	if (array1_size <= 1)
		return 1.0;

	if (array2_size <= 1)
		return 1.0;

	double fmean1 = 0.0, fmean2 = 0.0;
	fmean1 = avg (array1, array1_size);
	fmean2 = avg (array2, array2_size);

	if (fmean1 == fmean2)
		return 1.0;

	double usv1 = 0.0, usv2 = 0.0;

	for (size_t x = 0; x < array1_size; x++)
	{
		usv1 += (array1[x] - fmean1) * (array1[x] - fmean1);
	}

	for (size_t x = 0; x < array2_size; x++)
	{
		usv2 += (array2[x] - fmean2) * (array2[x] - fmean2);
	}

	usv1 = usv1 / (array1_size - 1);
	usv2 = usv2 / (array2_size - 1);
	const double welch_t_statistic = (fmean1 - fmean2) / sqrt (usv1 / array1_size + usv2 / array2_size);
	const double dof = pow ((usv1 / array1_size + usv2 / array2_size), 2.0) / ((usv1 * usv1) / (array1_size * array1_size * (array1_size - 1)) + (usv2 * usv2) / (array2_size * array2_size * (array2_size - 1)));

	return PfromT (welch_t_statistic, dof);
}

double PValuePaired (double *array1, double *array2, const size_t array_size)
{
	double *ABdiff = NULL;
	double mean;
	double std;
	double welch_t_statistic;
	double dof = array_size - 1;
	int i;
	ABdiff = malloc (array_size * sizeof (double));

	if (ABdiff == NULL)
		exit (EXIT_FAILURE);

	for (i = 0; i < array_size; i++)
		ABdiff[i] = array1[i] - array2[i];

	mean = avg (ABdiff, array_size);
	std = SDPop (ABdiff, array_size);
	welch_t_statistic = mean / (std / sqrt (array_size - 1));
	free (ABdiff);

	return PfromT (welch_t_statistic, dof);
}

double Lgamma(double x)
{
		double coef[6] = { 76.18009172947146,
				-86.50532032941677, 24.01409824083091,
				-1.231739572450155, 0.1208650973866179E-2,
				-0.5395239384953E-5 };
		double LogSqrtTwoPi = 0.91893853320467274178;
		double denom = x + 1;
		double y = x + 5.5;
		double series = 1.000000000190015;
		for (int i = 0; i < 6; ++i)
		{
				series += coef[i] / denom;
				denom += 1.0;
		}
		return (LogSqrtTwoPi + (x + 0.5) * log(y) - y + log(series / x));
}

struct pair getMinMax (double arr[], int n)
{
	struct pair minmax;
	int i;

	if (n % 2 == 0)
	{
		if (arr[0] > arr[1])
		{
			minmax.max = arr[0];
			minmax.min = arr[1];
		}
		else
		{
			minmax.min = arr[0];
			minmax.max = arr[1];
		}

		i = 2;
	}
	else
	{
		minmax.min = arr[0];
		minmax.max = arr[0];
		i = 1;
	}

	while (i < n - 1)
	{
		if (arr[i] > arr[i + 1])
		{
			if (arr[i] > minmax.max)
				minmax.max = arr[i];

			if (arr[i + 1] < minmax.min)
				minmax.min = arr[i + 1];
		}
		else
		{
			if (arr[i + 1] > minmax.max)
				minmax.max = arr[i + 1];

			if (arr[i] < minmax.min)
				minmax.min = arr[i];
		}

		i += 2;
	}

	return minmax;
}

static int getstr(char **lineptr, size_t *n, FILE *stream, char terminator, size_t offset)
{
	int nchars_avail;
	char *read_pos;
	int ret;

	if (!lineptr || !n || !stream)
		return -1;

	if (!*lineptr) 
	{
		*n = MIN_CHUNK;
		*lineptr = malloc(*n);
		if (!*lineptr)
			return -1;
	}

	nchars_avail = *n - offset;
	read_pos = *lineptr + offset;

	for (;;) 
	{
		int c = getc(stream);

		assert(*n - nchars_avail == read_pos - *lineptr);
		if (nchars_avail < 2) 
		{
			if (*n > MIN_CHUNK)
				*n *= 2;
			else
				*n += MIN_CHUNK;

			nchars_avail = *n + *lineptr - read_pos;
			*lineptr = realloc(*lineptr, *n);
		
			if (!*lineptr)
				return -1;
		
			read_pos = *n - nchars_avail + *lineptr;
			assert(*n - nchars_avail == read_pos - *lineptr);
		}

		if (c == EOF || ferror (stream)) 
		{
			if (read_pos == *lineptr)
				return -1;
			else
				break;
		}

		*read_pos++ = c;
		nchars_avail--;

		if (c == terminator)
			break;
	}

	*read_pos = '\0';

	ret = read_pos - (*lineptr + offset);
	return ret;
}

int readline (char **lineptr, size_t *n, FILE *stream)
{
	return getstr(lineptr, n, stream, '\n', 0);
}

int parsecmdline(struct cmdparams *cmdparams,int c, char*v[])
{
	int x;
	int f=0;
	int m=0;

	cmdparams->paired = 0;
	cmdparams->filter = 0;
	cmdparams->pmean = 0.0;
	cmdparams->clevel = .05;
	cmdparams->fname1 = NULL;
	cmdparams->fname2 = NULL;

	if (c < 2)
	{
		printf("\n\n./csvtop <parameters> <file1> [<file2>]\n");
		printf("\nParameters:");
		printf("\n    -p          Paired Data (Must have equal number of records in each file)");
		printf("\n    -c <level>  Confidence Level (0-100)");
		printf("\n    -f          Filter data using Chauvenet Outlier Removal before running test");
		printf("\n    -m <mean>   Predicted mean value for one sample test");
		
		return -1;
	}
			
	for (x=1;x<c;x++)
	{
		if (strstr(v[x],"-c")>0) // Confidence Level
		{
			cmdparams->clevel = (100.0 - atof(v[x+1])) / 100;

			if (cmdparams->clevel <= 0.0 || cmdparams->clevel >= 1.0)
			{
				printf ("\n\nConfidence Level Out of Bounds (0.0 < C < 100.0)\n\n");
				return -1;
			}
		}

		if (strstr(v[x],"-p")>0) // Paired
			cmdparams->paired=1;

		if (strstr(v[x],"-f")>0) // Chauvenet filter
			cmdparams->filter=1;

		if (strstr(v[x],"-m")>0) // Predicted Mean
		{
			m=1;
			cmdparams->pmean = atof(v[x+1]);
		}

		if (strstr(v[x],".csv")>0) // Filenames 
		{
			f++;
			switch (f)
			{
				case 1:
				{
					cmdparams->fname1 = strdup(v[x]);
					break;
				}

				case 2:	
				{
					cmdparams->fname2 = strdup(v[x]);
					break;
				}
			}
		}
	}

	if (cmdparams->paired == 1 && f!=2) 
	{
		printf ("\n\nPaired data specified but wrong number of files given!\n\n");
		return -1;
	}

	if (m == 1 && f!=1) 
	{
		printf ("\n\nPredicted Mean value specified but wrong number of files given!\n\n");
		return -1;
	}

	if (f==1 && m==0)
	{
		printf ("\n\nNo Predicted Mean value given!\n\n");
		return -1;
	}

	return 0;
}


// Calculate One sided P-Value
int OneSample(char* fname1, double clevel, double mean, int filter)
{
	char *line = NULL;
	FILE *in;
	int i, err = -1;;
	size_t len = 0;
	ssize_t read;
	double *bufferA = NULL;
	double *bufferAO = NULL;
	int countA = 0, lastCountA = 0;
	int numberOfFields;
	char **parsed = NULL;
	double p, avgA, SDA, SEA, SDAP;
	double sig2P, sig1P;
	double ccfs = 0.5;

	struct pair minmaxA;

	double Z;
	double cuA, clA;

	if ((in=fopen (fname1, "r")) == NULL)
		goto fail;

	while ((read = readline (&line, &len, in)) != -1)
	{
		if (!(parsed = CSVParse (line, &numberOfFields)))
		{
			printf ("String parsing failed!\n");
			goto fail;
		}

		countA += numberOfFields;

		bufferA = realloc (bufferA, sizeof (double) * countA);

		if (bufferA == NULL)
			goto fail;

		for (i = lastCountA; i < countA; i++)
		{
			bufferA[i] = atof (parsed[i - lastCountA]);
		} 

  	cleanupStrings (parsed, numberOfFields);

		lastCountA = countA;
	}

	fclose (in);
	in = NULL;

	printf ("\n\n\n%sP-Value criteria for FALSE null hypothesis < %s%.6g", KCYN, KYEL, clevel);
	printf ("\n\n%s *** Performing One Sample Test ***", KCYN);

	if (filter == 1)
	{
		printf ("\n\n\n\n%s *** Data after Chauvenets Criterion Outlier Removal Filter ***", KRED);	
		printf ("\n\n%sRemoving Outliers", KMAG);

		countA = RemoveOutliersUnpaired (bufferA, &bufferAO, lastCountA, ccfs);
	}
	else
	{
		printf ("\n\n%s *** Raw Data ***", KRED);
	}

	Z = critz(clevel / 2);
	avgA = avg(bufferA, countA);

	minmaxA = getMinMax(bufferA, countA);
	SDA = SDSamp(bufferA, countA);
	SDAP = SDPop(bufferA, countA);
	cuA = avgA + Z * (SDA / sqrt(countA));
	clA = avgA - Z * (SDA / sqrt(countA));
	SEA = StandardError(bufferA, countA);

	printf ("\n\n%sCount = %s%d", KGRN, KYEL, countA);
	printf ("\n\n%sMin = %s%.6g", KGRN, KYEL, minmaxA.min);
	printf ("\n%sMax = %s%.6g", KGRN, KYEL, minmaxA.max);

	printf ("\n\n%sHypothesis Mean = %s%.6g", KGRN, KYEL, mean);
	printf ("\n\n%sSample Mean A = %s%.6g", KGRN, KYEL, avgA);
	printf ("\n%sSample Median A = %s%.6g", KGRN, KYEL, median(bufferA, countA));

	if (mean < avgA)
	{
		printf ("\n\n%sSample A Mean Difference = %s+%.6g", KGRN, KYEL, fabs(avgA - mean));
		printf ("\n%sSample A Mean %% Change = %s+%.1f%%", KGRN, KYEL, fabs(PerDiff(mean, avgA)));
	}

	if (mean > avgA)
	{
		printf ("\n\n%sSample A Mean Difference = %s-%.6g", KGRN, KYEL, fabs(avgA - mean));
		printf ("\n%sSample A Mean %% Change = %s-%.1f%%", KGRN, KYEL, fabs(PerDiff(mean, avgA)));
	}

	printf ("\n\n%sA %0.0f%% CI = %s%.6g to %.6g", KGRN, (1.0 - clevel) * 100, KYEL, clA, cuA);
	
	printf ("\n\n%sSample SD A = %s%.6g", KGRN, KYEL, SDA);
	printf ("\n%sPopulation SD A = %s%.6g", KGRN, KYEL, SDAP);

	printf ("\n\n%sSample SE A = %s%.6g", KGRN, KYEL, SEA);


	printf("\n\n\n%s*** Welch t-test ***", KBLU);

	p = PValue(bufferA, countA, mean);
	sig2P = critz(p);
	sig1P = critz(p * 0.5);

	printf("\n%sNull Hypothesis is",KGRN);

	if (p <= clevel)
		printf("%s FALSE ",KCYN);
	else
		printf("%s TRUE ",KCYN);

	printf ("%sfor Two Sided test", KGRN);

	printf ("\n\n%sP-Value Two Sided = %s%.6g", KGRN, KYEL, p);
	printf ("\n%sSigma Level %s %s%1.1f", KGRN,(sig2P < 5.99) ? "=" : ">", KYEL, sig2P);

	printf("\n\n%sNull Hypothesis is",KGRN);

	if (0.5 * p <= clevel)
		printf("%s FALSE ",KCYN);
	else
		printf("%s TRUE ",KCYN);

	printf ("%sfor One Sided test", KGRN);

	if (avgA < mean)
		printf ("\n\n%sP-Value One Sided A < MEAN = %s%.6g", KGRN, KYEL, 0.5 * p);
	else
		printf ("\n\n%sP-Value One Sided A > MEAN = %s%.6g", KGRN, KYEL, 0.5 * p);
	printf ("\n%sSigma Level %s %s%1.1f", KGRN,(sig2P < 5.99) ? "=" : ">", KYEL, sig1P);

	err = 0;
	fail:

	free (line);
	free (bufferAO);
	free (bufferA);
	if (in) fclose (in);

	return err;
}

// Calculate Two sided P-Value
int TwoSample(char* fname1, char* fname2, double clevel, int filter, int paired)
{
	double p, avgA, avgB, SDA, SDB, SEA, SEB, SDAP, SDBP;
	double SED, SDD, cuAD, clAD, MoD;
	double sig2P, sig1P;
	double r, pr, tr, sr;
	double ccfs = 0.5;
	char *line = NULL;
	FILE *in;
	int i, err = -1;;
	size_t len = 0;
	ssize_t read;
	double *bufferA = NULL;
	double *bufferAO = NULL;
	double *bufferB = NULL;
	double *bufferBO = NULL;
	int countA = 0, lastCountA = 0;
	int countB = 0, lastCountB = 0;
	int numberOfFields;
	char **parsed = NULL;
	double Z;
	double cuA, clA, cuB, clB;

	struct pair minmaxA;
	struct pair minmaxB;

	/* Load A data */
	if ((in=fopen (fname1, "r")) == NULL)
		goto fail;

	while ((read = readline (&line, &len, in)) != -1)
	{
		if (!(parsed = CSVParse (line, &numberOfFields)))
		{
			printf ("String parsing failed!\n");
			goto fail;
		}

		countA += numberOfFields;

		bufferA = realloc (bufferA, sizeof (double) * countA);

		if (bufferA == NULL)
			goto fail;

		for (i = lastCountA; i < countA; i++)
		{
			bufferA[i] = atof (parsed[i - lastCountA]);
		} 

		cleanupStrings (parsed, numberOfFields);

		lastCountA = countA;
	}

	fclose (in);
	in = NULL;

	/* Load B data */
	if ((in=fopen (fname2, "r")) == NULL)
		goto fail;

	while ((read = readline (&line, &len, in)) != -1)
	{
		if (!(parsed = CSVParse (line, &numberOfFields)))
		{
			printf ("String parsing failed!\n");
			goto fail;
		}

		countB += numberOfFields;

		bufferB = realloc (bufferB, sizeof (double) * countB);

		if (bufferB == NULL)
			goto fail;

		for (i = lastCountB; i < countB; i++)
		{
			bufferB[i] = atof (parsed[i - lastCountB]);
		} 

		cleanupStrings (parsed, numberOfFields);

		lastCountB = countB;
	}

	fclose (in);
	in = NULL;

	if (paired == 1)
	{
		if (countA != countB)
		{
			printf("\n\nPaired data specified, but unequal number of A/B values used");
			goto fail;
		}
	}

	printf ("\n\n\n%sP-Value criteria for FALSE null hypothesis < %s%.6g", KCYN, KYEL, clevel);
	printf ("\n\n%s *** Performing Two Sample Test ***", KCYN);

	if (filter == 1)
	{
		printf ("\n\n\n\n%s *** Data after Chauvenets Criterion Outlier Removal Filter ***", KRED);	

		if (paired == 1)
		{
			printf ("\n\n%sRemoving Outliers in Pairs", KMAG);
			countA = RemoveOutliersPaired(bufferA, &bufferAO, bufferB, &bufferBO, countA, ccfs);
			countB = countA;
		}
		else
		{
			printf ("\n\n%sRemoving Outliers", KMAG);
			countA = RemoveOutliersUnpaired (bufferA, &bufferAO, lastCountA, ccfs);
			printf ("\n\n%sRemoving Outliers", KMAG);
			countB = RemoveOutliersUnpaired (bufferB, &bufferBO, lastCountB, ccfs);
		}
	}
	else
	{
		printf ("\n\n%s *** Raw Data ***", KRED);
	}
	
	Z = critz(clevel / 2);
	avgA = avg(bufferA, countA);
	avgB = avg(bufferB, countB);
	minmaxA = getMinMax(bufferA, countA);
	minmaxB = getMinMax(bufferB, countB);
	SDA = SDSamp(bufferA, countA);
	SDB = SDSamp(bufferB, countB);
	SDAP = SDPop(bufferA, countA);
	SDBP = SDPop(bufferB, countB);

	cuA = avgA + Z * (SDA / sqrt(countA));
	clA = avgA - Z * (SDA / sqrt(countA));
	cuB = avgB + Z * (SDB / sqrt(countB));
	clB = avgB - Z * (SDB / sqrt(countB));
	SEA = StandardError(bufferA, countA);
	SEB = StandardError(bufferB, countB);

	printf ("\n\n%sA Count = %s%d", KGRN, KYEL, countA);
	printf ("\n%sB Count = %s%d", KGRN, KYEL, countB);

	printf ("\n\n%sA Min = %s%.6g", KGRN, KYEL, minmaxA.min);
	printf ("\n%sA Max = %s%.6g", KGRN, KYEL, minmaxA.max);
	printf ("\n\n%sB Min = %s%.6g", KGRN, KYEL, minmaxB.min);
	printf ("\n%sB Max = %s%.6g", KGRN, KYEL, minmaxB.max);

	printf ("\n\n%sSample Mean A = %s%.6g", KGRN, KYEL, avgA);
	printf ("\n%sSample Mean B = %s%.6g", KGRN, KYEL, avgB);
	printf ("\n\n%sSample Median A = %s%.6g", KGRN, KYEL, median(bufferA, countA));
	printf ("\n%sSample Median B = %s%.6g", KGRN, KYEL, median(bufferB, countB));

	if (avgA < avgB)
	{
		printf ("\n\n%sSample Mean Difference = %s+%.6g", KGRN, KYEL, fabs(avgB - avgA));
		printf ("\n%sSample Mean %% Change = %s+%.1f%%", KGRN, KYEL, fabs(PerDiff(avgA, avgB)));
	}

	if (avgA > avgB)
	{
		printf ("\n\n%sSample Mean Difference = %s-%.6g", KGRN, KYEL, fabs(avgB - avgA));
		printf ("\n%sSample Mean %% Change = %s-%.1f%%", KGRN, KYEL, fabs(PerDiff(avgA, avgB)));
	}

	printf ("\n\n%sA %0.0f%% CI = %s%.6g to %.6g", KGRN, (1.0 - clevel) * 100, KYEL, clA, cuA);
	printf ("\n%sB %0.0f%% CI = %s%.6g to %.6g", KGRN, (1.0 - clevel) * 100, KYEL, clB, cuB);

	printf ("\n\n%sSample SD A = %s%.6g", KGRN, KYEL, SDA);
	printf ("\n%sPopulation SD A = %s%.6g", KGRN, KYEL, SDAP);

	printf ("\n\n%sSample SD B = %s%.6g", KGRN, KYEL, SDB);
	printf ("\n%sPopulation SD B = %s%.6g", KGRN, KYEL, SDBP);

	if (SDA < SDB)
	{
		printf("\n\n%sSample SD Difference = %s+%.6g", KGRN, KYEL,fabs(SDB - SDA));
		printf("\n%sSample SD %% Change = %s+%.1f%%", KGRN, KYEL,fabs(PerDiff(SDA, SDB)));
	}

	if (SDA > SDB)
	{
		printf("\n\n%sSample SD Difference = %s-%.6g", KGRN, KYEL,fabs(SDB - SDA));
		printf("\n%sSample SD %% Change = %s-%.1f%%", KGRN, KYEL,fabs(PerDiff(SDA, SDB)));
	}

	if (paired == 1)
	{
		if (SEofDifferences(&SED, bufferA, bufferB, countA)<0) goto fail;
		if (SDofDifferences(&SDD, bufferA, bufferB, countA)<0) goto fail;
		if (MeanofDifferences(&MoD, bufferA, bufferB, countA)<0) goto fail;

		cuAD = MoD + Z * (SDD / sqrt(countA));
		clAD = MoD - Z * (SDD / sqrt(countA));

		printf("\n\n%sSD of Sample Differences = %s%.6g", KGRN, KYEL, SDD);
		printf("\n%sSE of Sample Differences = %s%.6g", KGRN, KYEL, SED);
		printf ("\n%sSample Differences %0.0f%% CI = %s%.6g to %.6g", KGRN, (1.0 - clevel) * 100, KYEL, clAD, cuAD);
	}

	printf ("\n\n%sSample SE A = %s%.6g", KGRN, KYEL, SEA);
	printf ("\n%sSample SE B = %s%.6g", KGRN, KYEL, SEB);

	if (paired == 0)
	{
		printf("\n\n\n%s*** Welch t-test UnPaired***", KBLU);

		p = PValueUnpaired(bufferA, countA, bufferB, countB);
		sig2P = critz(p);
		sig1P = critz(p * 0.5);

		printf("\n%sNull Hypothesis is",KGRN);

		if (p <= clevel)
			printf("%s FALSE ",KCYN);
		else
			printf("%s TRUE ",KCYN);

		printf ("%sfor Two Sided test", KGRN);

		printf ("\n\n%sP-Value Two Sided = %s%.6g", KGRN, KYEL, p);
		printf ("\n%sSigma Level %s %s%1.1f", KGRN,(sig2P < 5.99) ? "=" : ">", KYEL, sig2P);

		printf("\n\n%sNull Hypothesis is",KGRN);

		if (0.5 * p <= clevel)
			printf("%s FALSE ",KCYN);
		else
			printf("%s TRUE ",KCYN);

		printf ("%sfor One Sided test", KGRN);
		if (avgA < avgB)
			printf ("\n\n%sP-Value One Sided A < B = %s%.6g", KGRN, KYEL, 0.5 * p);
		else
			printf ("\n\n%sP-Value One Sided A > B = %s%.6g", KGRN, KYEL, 0.5 * p);
		printf ("\n%sSigma Level %s %s%1.1f", KGRN,(sig2P < 5.99) ? "=" : ">", KYEL, sig1P);
	}
	else
	{
		printf("\n\n\n%s*** Welch t-test Paired***", KBLU);

		p = PValuePaired(bufferA, bufferB, countA);
		sig2P = critz(p);
		sig1P = critz(p * 0.5);

		printf("\n%sNull Hypothesis is",KGRN);

		if (p <= clevel)
			printf("%s FALSE ",KCYN);
		else
			printf("%s TRUE ",KCYN);

		printf ("%sfor Two Sided test", KGRN);

		printf ("\n\n%sP-Value Two Sided = %s%.6g", KGRN, KYEL, p);
		printf ("\n%sSigma Level %s %s%1.1f", KGRN,(sig2P < 5.99) ? "=" : ">", KYEL, sig2P);

		printf("\n\n%sNull Hypothesis is",KGRN);

		if (0.5 * p <= clevel)
			printf("%s FALSE ",KCYN);
		else
			printf("%s TRUE ",KCYN);

		printf ("%sfor One Sided test", KGRN);
		if (avgA < avgB)
			printf ("\n\n%sP-Value One Sided A < B = %s%.6g", KGRN, KYEL, 0.5 * p);
		else
			printf ("\n\n%sP-Value One Sided A > B = %s%.6g", KGRN, KYEL, 0.5 * p);
		printf ("\n%sSigma Level %s %s%1.1f", KGRN,(sig2P < 5.99) ? "=" : ">", KYEL, sig1P);
	}

	printf("\n\n\n%s*** Pearson Correlation Coefficient ***", KBLU);

	r = R(bufferA, bufferB, countA);
	tr = r / sqrt((1 - r * r) / (countA - 2));
	pr = PfromT(tr, countA - 2);

	printf ("\n%sA to B has ", KGRN);

	if (r == 0.0) printf ("%sNo", KCYN);
	if (r == 1.0) printf ("%sPerfect Positive", KCYN);
	if (r == -1.0) printf ("%sPerfect Negative", KCYN);

	if (r > 0.0 && r < 0.3) printf ("%sWeak Positive", KCYN);
	if (r >= 0.3 && r < 0.7) printf ("%sModerate Positive", KCYN);
	if (r >= 0.7 && r < 1.00) printf ("%sStrong Positive", KCYN);

	if (r < 0.0 && r > -0.3) printf ("%sWeak Negative", KCYN);
	if (r <= -0.3 && r > -0.7) printf ("%sWeak Negative", KCYN);
	if (r <= -0.7 && r > -1.00) printf ("%sWeak Negative", KCYN);

	printf ("%s Correlation", KGRN);

	printf ("\n\n%sR-Value = %s%.6g", KGRN, KYEL, r);
	printf ("\n\n%sCoefficient of Determination = %s%.6g", KGRN, KYEL, r * r);
	printf ("\n\n%sP-Value = %s%.6g", KGRN, KYEL, pr);

	sr = critz(pr);
	printf ("\n%sSigma Level %s %s%1.1f", KGRN,(sig2P < 5.99) ? "=" : ">", KYEL, sr);

	if (pr <= clevel)
		printf("\n\n%sP-Value is %sSignificant",KGRN,KCYN);
	else
		printf("\n\n%sP-Value is %sNot Significant ",KGRN,KCYN);

	err = 0;
	fail:

	free (line);
	free (bufferA);
	free (bufferAO);
	free (bufferB);
	free (bufferBO);

	if (in) fclose (in);
	return err;
}

int main (int argc, char *argv[])
{
	int err=-1;
	struct cmdparams cmdparams;		

	if (parsecmdline(&cmdparams, argc, argv)<0) goto done; // Parse the command line parameters

	if (cmdparams.fname2 == NULL) // Choose One or Two samples
		OneSample(cmdparams.fname1, cmdparams.clevel, cmdparams.pmean, cmdparams.filter);
	else
		TwoSample(cmdparams.fname1, cmdparams.fname2, cmdparams.clevel, cmdparams.filter, cmdparams.paired);
	
	err=0;

	done:

	/* Clean up after ourselves */
	printf ("%s\n\n", KNRM);
	
	free(cmdparams.fname1);
	free(cmdparams.fname2);

	printf ("\n\n");

	return err;
}
