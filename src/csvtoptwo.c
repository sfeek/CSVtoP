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

int RemoveOutliers (double *in, double **out, int n, double sensitivity)
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

double PfromT (const double welch_t_statistic, const double dof)
{
	const double a = dof / 2;
	double value = dof / (welch_t_statistic * welch_t_statistic + dof);

	if ((isinf (value) != 0) || (isnan (value) != 0))
		return 1.0;

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

int main (int argc, char *argv[])
{
	char **parsed = NULL;
	FILE *in;
	int numberOfFields;
	char *line = NULL;
	int i;
	size_t len = 0;
	ssize_t read;
	double *bufferA = NULL;
	double *bufferB = NULL;
	double *bufferAO = NULL;
	double *bufferBO = NULL;
	int countA = 0, lastCountA = 0;
	int countB = 0, lastCountB = 0;
	double p, avgA, avgB, SDA, SDB, clevel;
	double sig2P, sig1P;
	struct pair minmaxA;
	struct pair minmaxB;
	double Z;
	double cuA, clA, cuB, clB;

	if (argc != 4)
	{
		printf ("\n\n Usage: csvtoptwo <confidence_level> <fname A>.csv <fname B>.csv\n");
		exit (EXIT_FAILURE);
	}

	/* get Confidence level */
	clevel = atof (argv[1]);

	if (clevel <= 0.0 || clevel >= 1.0)
	{
		printf ("\n\nConfidence Level Out of Bounds (0.0 < C < 1.0)\n\n");
		exit (EXIT_FAILURE);
	}

	/* Read in A array */
	/* Open the files */
	in = fopen (argv[2], "r");

	if (in == NULL)
		exit (EXIT_FAILURE);


	/* Read each line */
	while ((read = readline (&line, &len, in)) != -1)
	{
		/* Parse it! */
		if (! (parsed = CSVParse (line, &numberOfFields)))
		{
			printf ("String parsing failed!\n");
			return 1;
		}

		/* Keep track of how many field total */
		countA += numberOfFields;
		/* Make Space */
		bufferA = realloc (bufferA, sizeof (double) * countA);

		if (bufferA == NULL)
			exit (EXIT_FAILURE);

		/* Fill it! */
		for (i = lastCountA; i < countA; i++)
		{
			bufferA[i] = atof (parsed[i - lastCountA]);
		}

		/* Cleanup parse fields */
		cleanupStrings (parsed, numberOfFields);
		/* Keep track of our last count */
		lastCountA = countA;
	}

	/* Clean up and close out */
	free (line);
	line = NULL;
	fclose (in);
	/* Read in B array */
	/* Open the files */
	in = fopen (argv[3], "r");

	if (in == NULL)
		exit (EXIT_FAILURE);

	/* Read each line */
	while ((read = readline (&line, &len, in)) != -1)
	{
		/* Parse it! */
		if (! (parsed = CSVParse (line, &numberOfFields)))
		{
			printf ("String parsing failed!\n");
			return 1;
		}

		/* Keep track of how many field total */
		countB += numberOfFields;
		/* Make Space */
		bufferB = realloc (bufferB, sizeof (double) * countB);

		if (bufferB == NULL)
			exit (EXIT_FAILURE);

		/* Fill it! */
		for (i = lastCountB; i < countB; i++)
		{
			bufferB[i] = atof (parsed[i - lastCountB]);
		}

		/* Cleanup parse fields */
		cleanupStrings (parsed, numberOfFields);
		/* Keep track of our last count */
		lastCountB = countB;
	}

	/* Clean up and close out */
	free (line);
	line = NULL;
	fclose (in);
	/* Show the results */
	Z = critz (clevel / 2);
	avgA = avg (bufferA, lastCountA);
	avgB = avg (bufferB, lastCountB);
	minmaxA = getMinMax (bufferA, lastCountA);
	minmaxB = getMinMax (bufferB, lastCountB);
	SDA = SDSamp (bufferA, lastCountA);
	SDB = SDSamp (bufferB, lastCountB);
	cuA = avgA + Z * (SDA / sqrt (countA));
	clA = avgA - Z * (SDA / sqrt (countA));
	cuB = avgB + Z * (SDB / sqrt (countB));
	clB = avgB - Z * (SDB / sqrt (countB));
	printf ("\n\n\n%sP-Value criteria for FALSE null hypothesis < %s%.6g", KCYN, KYEL, clevel);
	printf ("\n\n%s *** Raw Data ***", KRED);
	printf ("\n\n%sA Count = %s%d", KGRN, KYEL, lastCountA);
	printf ("\n%sB Count = %s%d\n", KGRN, KYEL, lastCountB);
	printf ("\n%sA Min = %s%.6g", KGRN, KYEL, minmaxA.min);
	printf ("\n%sA Max = %s%.6g", KGRN, KYEL, minmaxA.max);
	printf ("\n%sB Min = %s%.6g", KGRN, KYEL, minmaxB.min);
	printf ("\n%sB Max = %s%.6g\n", KGRN, KYEL, minmaxB.max);
	printf ("\n%sSample Mean A = %s%.6g", KGRN, KYEL, avgA);
	printf ("\n%sSample Mean B = %s%.6g", KGRN, KYEL, avgB);

	if (avgA < avgB)
	{
		printf ("\n%sSample Mean Difference = %s+%.6g", KGRN, KYEL, fabs (avgB - avgA));
		printf ("\n%sSample Mean %% Change = %s+%.2g%%", KGRN, KYEL, PerDiff (avgA, avgB));
	}

	if (avgA > avgB)
	{
		printf ("\n%sSample Mean Difference = %s-%.6g", KGRN, KYEL, fabs (avgB - avgA));
		printf ("\n%sSample Mean %% Change = %s-%.2g%%", KGRN, KYEL, PerDiff (avgA, avgB));
	}

	printf ("\n\n%sA %.2g%% CI = %s %.6g - %.6g", KGRN, (1.0 - clevel) * 100, KYEL, clA, cuA);
	printf ("\n%sB %.2g%% CI = %s %.6g - %.6g", KGRN, (1.0 - clevel) * 100, KYEL, clB, cuB);
	printf ("\n\n%sSample SD A = %s%.6g", KGRN, KYEL, SDA);
	printf ("\n%sSample SD B = %s%.6g", KGRN, KYEL, SDB);

	if (SDA < SDB)
	{
		printf ("\n%sSample SD Difference = %s+%.6g", KGRN, KYEL, fabs (SDB - SDA));
		printf ("\n%sSample SD %% Change = %s+%.2g%%", KGRN, KYEL, PerDiff (SDA, SDB));
	}

	if (SDA > SDB)
	{
		printf ("\n%sSample SD Difference = %s-%.6g", KGRN, KYEL, fabs (SDB - SDA));
		printf ("\n%sSample SD %% Change = %s-%.2g%%", KGRN, KYEL, PerDiff (SDA, SDB));
	}

	printf ("\n\n\n%s*** Welch t-test Unpaired ***", KBLU);
	p = PValueUnpaired (bufferA, lastCountA, bufferB, lastCountB);
	sig2P = critz (p);
	sig1P = critz (p * 0.5);

	if (p <= clevel)
		printf ("\n%sNull Hypothesis is %sFALSE%s for Two Sided test", KGRN, KCYN, KGRN);
	else
		printf ("\n%sNull Hypothesis is %sTRUE%s for Two Sided test", KGRN, KCYN, KGRN);

	printf ("\n%sP-Value Two Sided = %s%.6g", KGRN, KYEL, p);
	printf ("\n%sSigma Level %s %s%3.1f \n", KGRN, (sig2P < 5.99) ? "=" : ">", KYEL, sig2P);

	if (avgA < avgB)
	{
		if (0.5 * p <= clevel)
			printf ("\n%sNull Hypothesis is %sFALSE%s for One Sided test", KGRN, KCYN, KGRN);
		else
			printf ("\n%sNull Hypothesis is %sTRUE%s for One Sided test", KGRN, KCYN, KGRN);

		printf ("\n%sP-Value One Sided A < B = %s%.6g", KGRN, KYEL, 0.5 * p);
		printf ("\n%sSigma Level %s %s%3.1f \n", KGRN, (sig1P < 5.99) ? "=" : ">", KYEL, sig1P);
	}
	else
	{
		if (0.5 * p <= clevel)
			printf ("\n%sNull Hypothesis is %sFALSE%s for One Sided test", KGRN, KCYN, KGRN);
		else
			printf ("\n%sNull Hypothesis is %sTRUE%s for One Sided test", KGRN, KCYN, KGRN);

		printf ("\n%sP-Value One Sided A > B = %s%.6g", KGRN, KYEL, 0.5 * p);
		printf ("\n%sSigma Level %s %s%3.1f \n", KGRN, (sig1P < 5.99) ? "=" : ">", KYEL, sig1P);
	}

	if (lastCountA == lastCountB)
	{
		printf ("\n\n\n%s*** Welch t-test Paired ***", KBLU);
		p = PValuePaired (bufferA, bufferB, lastCountA);
		sig2P = critz (p);
		sig1P = critz (p * 0.5);

		if (p <= clevel)
			printf ("\n%sNull Hypothesis is %sFALSE%s for Two Sided test", KGRN, KCYN, KGRN);
		else
			printf ("\n%sNull Hypothesis is %sTRUE%s for Two Sided test", KGRN, KCYN, KGRN);

		printf ("\n%sP-Value Two Sided = %s%.6g", KGRN, KYEL, p);
		printf ("\n%sSigma Level %s %s%3.1f \n", KGRN, (sig2P < 5.99) ? "=" : ">", KYEL, sig2P);

		if (avgA < avgB)
		{
			if (0.5 * p <= clevel)
				printf ("\n%sNull Hypothesis is %sFALSE%s for One Sided test", KGRN, KCYN, KGRN);
			else
				printf ("\n%sNull Hypothesis is %sTRUE%s for One Sided test", KGRN, KCYN, KGRN);

			printf ("\n%sP-Value One Sided A < B = %s%.6g", KGRN, KYEL, 0.5 * p);
			printf ("\n%sSigma Level %s %s%3.1f \n", KGRN, (sig1P < 5.99) ? "=" : ">", KYEL, sig1P);
		}
		else
		{
			if (0.5 * p <= clevel)
				printf ("\n%sNull Hypothesis is %sFALSE%s for One Sided test", KGRN, KCYN, KGRN);
			else
				printf ("\n%sNull Hypothesis is %sTRUE%s for One Sided test", KGRN, KCYN, KGRN);

			printf ("\n%sP-Value One Sided A > B = %s%.6g", KGRN, KYEL, 0.5 * p);
			printf ("\n%sSigma Level %s %s%3.1f \n", KGRN, (sig1P < 5.99) ? "=" : ">", KYEL, sig1P);
		}
	}

	printf ("\n\n\n\n%s *** Data after Chauvenets Criterion Outlier Removal Filter ***", KRED);
	printf ("\n\n%sRemoving Outliers from %sA", KMAG, KYEL);
	lastCountA = RemoveOutliers (bufferA, &bufferAO, lastCountA, 0.5);
	printf ("\n\n%sRemoving Outliers from %sB", KMAG, KYEL);
	lastCountB = RemoveOutliers (bufferB, &bufferBO, lastCountB, 0.5);
	/* Show the results */
	avgA = avg (bufferAO, lastCountA);
	avgB = avg (bufferBO, lastCountB);
	minmaxA = getMinMax (bufferAO, lastCountA);
	minmaxB = getMinMax (bufferBO, lastCountB);
	SDA = SDSamp (bufferAO, lastCountA);
	SDB = SDSamp (bufferBO, lastCountB);
	cuA = avgA + Z * (SDA / sqrt (countA));
	clA = avgA - Z * (SDA / sqrt (countA));
	cuB = avgB + Z * (SDB / sqrt (countB));
	clB = avgB - Z * (SDB / sqrt (countB));
	printf ("\n\n%sA Count = %s%d", KGRN, KYEL, lastCountA);
	printf ("\n%sB Count = %s%d\n", KGRN, KYEL, lastCountB);
	printf ("\n%sA Min = %s%.6g", KGRN, KYEL, minmaxA.min);
	printf ("\n%sA Max = %s%.6g", KGRN, KYEL, minmaxA.max);
	printf ("\n%sB Min = %s%.6g", KGRN, KYEL, minmaxB.min);
	printf ("\n%sB Max = %s%.6g\n", KGRN, KYEL, minmaxB.max);
	printf ("\n%sSample Mean A = %s%.6g", KGRN, KYEL, avgA);
	printf ("\n%sSample Mean B = %s%.6g", KGRN, KYEL, avgB);

	if (avgA < avgB)
	{
		printf ("\n%sSample Mean Difference = %s+%.6g", KGRN, KYEL, fabs (avgB - avgA));
		printf ("\n%sSample Mean %% Change = %s+%.2g%%", KGRN, KYEL, PerDiff (avgA, avgB));
	}

	if (avgA > avgB)
	{
		printf ("\n%sSample Mean Difference = %s-%.6g", KGRN, KYEL, fabs (avgB - avgA));
		printf ("\n%sSample Mean %% Change = %s-%.2g%%", KGRN, KYEL, PerDiff (avgA, avgB));
	}

	printf ("\n\n%sA %.2g%% CI = %s %.6g - %.6g", KGRN, (1.0 - clevel) * 100, KYEL, clA, cuA);
	printf ("\n%sB %.2g%% CI = %s %.6g - %.6g", KGRN, (1.0 - clevel) * 100, KYEL, clB, cuB);
	printf ("\n\n%sSample SD A = %s%.6g", KGRN, KYEL, SDA);
	printf ("\n%sSample SD B = %s%.6g", KGRN, KYEL, SDB);

	if (SDA < SDB)
	{
		printf ("\n%sSample SD Difference = %s+%.6g", KGRN, KYEL, fabs (SDB - SDA));
		printf ("\n%sSample SD %% Change = %s+%.2g%%", KGRN, KYEL, PerDiff (SDA, SDB));
	}

	if (SDA > SDB)
	{
		printf ("\n%sSample SD Difference = %s-%.6g", KGRN, KYEL, fabs (SDB - SDA));
		printf ("\n%sSample SD %% Change = %s-%.2g%%", KGRN, KYEL, PerDiff (SDA, SDB));
	}

	printf ("\n\n\n%s*** Welch t-test Unpaired ***", KBLU);
	p = PValueUnpaired (bufferAO, lastCountA, bufferBO, lastCountB);
	sig2P = critz (p);
	sig1P = critz (p * 0.5);

	if (p <= clevel)
		printf ("\n%sNull Hypothesis is %sFALSE%s for Two Sided test", KGRN, KCYN, KGRN);
	else
		printf ("\n%sNull Hypothesis is %sTRUE%s for Two Sided test", KGRN, KCYN, KGRN);

	printf ("\n%sP-Value Two Sided = %s%.6g", KGRN, KYEL, p);
	printf ("\n%sSigma Level %s %s%3.1f \n", KGRN, (sig2P < 5.99) ? "=" : ">", KYEL, sig2P);

	if (avgA < avgB)
	{
		if (0.5 * p <= clevel)
			printf ("\n%sNull Hypothesis is %sFALSE%s for One Sided test", KGRN, KCYN, KGRN);
		else
			printf ("\n%sNull Hypothesis is %sTRUE%s for One Sided test", KGRN, KCYN, KGRN);

		printf ("\n%sP-Value One Sided A < B = %s%.6g", KGRN, KYEL, 0.5 * p);
		printf ("\n%sSigma Level %s %s%3.1f \n", KGRN, (sig1P < 5.99) ? "=" : ">", KYEL, sig1P);
	}
	else
	{
		if (0.5 * p <= clevel)
			printf ("\n%sNull Hypothesis is %sFALSE%s for One Sided test", KGRN, KCYN, KGRN);
		else
			printf ("\n%sNull Hypothesis is %sTRUE%s for One Sided test", KGRN, KCYN, KGRN);

		printf ("\n%sP-Value One Sided A > B = %s%.6g", KGRN, KYEL, 0.5 * p);
		printf ("\n%sSigma Level %s %s%3.1f \n", KGRN, (sig1P < 5.99) ? "=" : ">", KYEL, sig1P);
	}

	/* Clean up after ourselves */
	printf ("%s\n\n", KNRM);
	free (bufferA);
	bufferA = NULL;
	free (bufferB);
	bufferB =  NULL;
	free (bufferAO);
	bufferAO = NULL;
	free (bufferBO);
	bufferBO = NULL;
	printf ("\n\n");
	/* And we are out of here! */
	return 0;
}
