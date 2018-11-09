#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CSVLib.h>
#include <math.h>

#define MAXLOG 7.09782712893383996732E2

double erf_local(double x);
double erfc_local(double a);

struct pair  
{ 
  double min; 
  double max; 
};   
  
double avg(double *buffer, int count)
{
	int i;
	double total=0;
	
	for (i=0;i<count;i++) total += buffer[i];

	return total/count;
}

double SDPop(double *buffer, int count)
{
	double mean, standardDeviation = 0.0;
	int i;

	mean = avg(buffer,count);

	for(i=0; i<count; ++i) standardDeviation += pow(buffer[i] - mean, 2);

	return sqrt(standardDeviation/count);
}

double SDSamp(double *buffer, int count)
{
	double mean, standardDeviation = 0.0;
	int i;

	mean = avg(buffer,count);

	for(i=0; i<count; ++i) standardDeviation += pow(buffer[i] - mean, 2);

	return sqrt(standardDeviation/(count-1));
}

int RemoveOutliers (double *in, double **out, int n, double sensitivity)
{
	int c = 0;
	int i;
	double test;
	double a = avg (in,n);
	double std = SDSamp (in,n);

	*out = malloc(n * sizeof (double));
	if (*out == NULL) return 0;

	for (i = 0;i < n; i++)
	{
		test = n * erfc(fabs(in[i] - a) / std);
		if (test < sensitivity) 
			printf ("\nThrew out #%d -> %.6g",i+1,in[i]);
		else
			(*out)[c++] = in[i];
	}
	
	return c;  
}

double PfromT (const double welch_t_statistic, const double dof)
{ 	
	const double a = dof/2;
	double value = dof/(welch_t_statistic*welch_t_statistic+dof);
	
	if ((isinf(value) != 0) || (isnan(value) != 0)) return 1.0;

	printf ("\nT = %.6g",welch_t_statistic);
	printf ("\ndf = %.6g\n",dof);

	const double beta = lgammal(a)+0.57236494292470009-lgammal(a+0.5);
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
 
	if (value < 0.0 || 1.0 < value) return value;
  	
	if (value == 0.0 || value == 1.0) return value;
  	
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
  
	if (ns == 0) rx = xx;
 
  	for (;;)
  	{
    		term = term * temp * rx / (pp + ai);
    		value = value + term;;
    		temp = fabs (term);
 
    		if (temp <= acu && temp <= acu * value)
    		{
      			value = value * exp (pp * log (xx) + (qq - 1.0) * log (cx) - beta) / pp;
 
      			if (indx) value = 1.0 - value;
      			break;
    		}
 
    		ai = ai + 1.0;
    		ns = ns - 1;
 
    		if (0 <= ns)
    		{
      			temp = qq - ai;
      			if (ns == 0) rx = xx;
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
	if (array1_size <= 1) return 1.0;
	if (array2_size <= 1) return 1.0;
	
	double fmean1 = 0.0, fmean2 = 0.0;
	
	fmean1 = avg(array1,array1_size);
	fmean2 = avg(array2,array2_size);
	
	if (fmean1 == fmean2) return 1.0;

	double usv1 = 0.0, usv2 = 0.0;

	for (size_t x = 0; x < array1_size; x++) 
	{
		usv1 += (array1[x]-fmean1)*(array1[x]-fmean1);
	}

	for (size_t x = 0; x < array2_size; x++) 
	{
		usv2 += (array2[x]-fmean2)*(array2[x]-fmean2);
	}
	
	usv1 = usv1/(array1_size-1);
	usv2 = usv2/(array2_size-1);

	const double welch_t_statistic = (fmean1-fmean2)/sqrt(usv1/array1_size+usv2/array2_size);
	const double dof = pow((usv1/array1_size+usv2/array2_size),2.0)/((usv1*usv1)/(array1_size*array1_size*(array1_size-1))+(usv2*usv2)/(array2_size*array2_size*(array2_size-1)));

	return PfromT(welch_t_statistic,dof);
}

double PValuePaired (double *array1, double *array2, const size_t array_size) 
{
	double *ABdiff = NULL;
	double mean;
	double std;
	double welch_t_statistic;
	double dof = array_size - 1;

	int i;
	
	ABdiff = malloc (array_size * sizeof(double)); 
	if (ABdiff == NULL) exit (EXIT_FAILURE);

	for (i=0;i<array_size;i++) ABdiff[i] = array1[i] - array2[i];

	mean = avg(ABdiff,array_size);

	std = SDPop(ABdiff,array_size);

	welch_t_statistic = mean/(std/sqrt(array_size-1));

	free (ABdiff);

	return PfromT(welch_t_statistic,dof);
}

struct pair getMinMax(double arr[], int n) 
{ 
	struct pair minmax;      
	int i;   
  
	if (n%2 == 0) 
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
    
	while (i < n-1)   
	{           
		if (arr[i] > arr[i+1])           
		{ 
			if(arr[i] > minmax.max)         
			minmax.max = arr[i]; 
			if(arr[i+1] < minmax.min)           
			minmax.min = arr[i+1];         
		}  
		else         
		{ 
			if (arr[i+1] > minmax.max)         
				minmax.max = arr[i+1]; 
			if (arr[i] < minmax.min)           
				minmax.min = arr[i];         
		}         
		i += 2;   
	}             
  
	return minmax; 
}   

int main (int argc, char *argv[])
{
	char **parsed=NULL;
	FILE *in;
	int numberOfFields;
	char *line=NULL;
	int i;
	size_t len = 0;
	ssize_t read;
	double *bufferA=NULL;
	double *bufferB=NULL;
	double *bufferAO=NULL;
	double *bufferBO=NULL;
	int countA=0,lastCountA=0;
	int countB=0,lastCountB=0;
	double p,avgA,avgB;
	struct pair minmaxA;
	struct pair minmaxB;

	if (argc != 3)
	{
		printf ("\n\n Usage: csvtop <fname A.csv> <fname B>.csv\n");
		exit (EXIT_FAILURE);
	}

	/* Read in A array */
	/* Open the files */
	in = fopen (argv[1], "r");

	if (in == NULL)
		exit (EXIT_FAILURE);

	/* Read each line */
	while ((read = getline (&line, &len, in)) != -1)
   	{
		/* Parse it! */
   		if (!(parsed = CSVParse (line,&numberOfFields)))
   		{
       		printf ("String parsing failed!\n");
       		return 1;
   		}

		/* Keep track of how many field total */
		countA += numberOfFields;

		/* Make Space */	
		bufferA = realloc (bufferA, sizeof(double) * countA);
		if (bufferA == NULL) exit (EXIT_FAILURE);

		/* Fill it! */
   		for (i=lastCountA;i<countA;i++)
   		{
       		bufferA[i] = atof (parsed[i-lastCountA]);
   		}

		/* Cleanup parse fields */
   	   	cleanupStrings (parsed,numberOfFields);

		/* Keep track of our last count */
		lastCountA = countA;
	}

	/* Clean up and close out */
	free (line);
	line=NULL;
	fclose (in);

	/* Read in B array */
	/* Open the files */
	in = fopen (argv[2], "r");

	if (in == NULL)
		exit (EXIT_FAILURE);

	/* Read each line */
	while ((read = getline (&line, &len, in)) != -1)
   	{
		/* Parse it! */
   		if (!(parsed = CSVParse (line,&numberOfFields)))
   		{
       		printf ("String parsing failed!\n");
       		return 1;
   		}

		/* Keep track of how many field total */
		countB += numberOfFields;

		/* Make Space */	
		bufferB = realloc (bufferB, sizeof(double) * countB);
		if (bufferB == NULL) exit (EXIT_FAILURE);

		/* Fill it! */
   		for (i=lastCountB;i<countB;i++)
   		{
       		bufferB[i] = atof (parsed[i-lastCountB]);
   		}

		/* Cleanup parse fields */
   	   	cleanupStrings (parsed,numberOfFields);

		/* Keep track of our last count */
		lastCountB = countB;
	}

	/* Clean up and close out */
	free (line);
	line=NULL;
	fclose (in);

	/* Show the results */
	avgA = avg (bufferA,lastCountA);
	avgB = avg (bufferB,lastCountB);
	minmaxA = getMinMax(bufferA,lastCountA);
	minmaxB = getMinMax(bufferB,lastCountB);
	
	printf ("\n\n *** Raw Data  ***");

	printf ("\n\nA Count = %d",lastCountA);
	printf ("\nB Count = %d\n",lastCountB);
	printf ("\nA Min = %.6g", minmaxA.min);
	printf ("\nA Max = %.6g", minmaxA.max);
	printf ("\nB Min = %.6g", minmaxB.min);
	printf ("\nB Max = %.6g\n", minmaxB.max);
	printf ("\nAVG A = %.6g ",avgA);
	printf ("\nAVG B = %.6g \n",avgB);
	printf ("\nSD A = %.6g ",SDSamp(bufferA,lastCountA));
	printf ("\nSD B = %.6g \n",SDSamp(bufferB,lastCountB));



	printf ("\n*** Welch t-test Unpaired ***");

	p = PValueUnpaired (bufferA,lastCountA,bufferB,lastCountB);

	printf ("\nP-Value Two Sided = %.6g ",p);

	if (avgA <= avgB)
	{
		printf ("\nP-Value One Sided A < B = %.6g ", 0.5*p);
		printf ("\nP-Value One Sided A > B = %.6g ", 1 - 0.5*p);
	}
	else
	{
		printf ("\nP-Value One Sided A < B = %.6g ", 1 - 0.5*p);
		printf ("\nP-Value One Sided A > B = %.6g ", 0.5*p);
	}

	if (lastCountA == lastCountB)
	{
		printf ("\n\n*** Welch t-test Paired ***");
		
		p = PValuePaired(bufferA,bufferB,lastCountA);
	
		printf ("\nP-Value Two Sided = %.6g ",p);

		if (avgA <= avgB)
		{
			printf ("\nP-Value One Sided A < B = %.6g ", 0.5*p);
			printf ("\nP-Value One Sided A > B = %.6g ", 1 - 0.5*p);
		}
		else
		{
			printf ("\nP-Value One Sided A < B = %.6g ", 1 - 0.5*p);
			printf ("\nP-Value One Sided A > B = %.6g ", 0.5*p);
		}
	}
	
	printf ("\n\n\n *** Chauvenets Criterion Outlier Removal ***");

	printf ("\n\nRemoving Outliers from A");
	lastCountA = RemoveOutliers (bufferA,&bufferAO,lastCountA,0.5);
	
	printf ("\n\nRemoving Outliers from B");
	lastCountB = RemoveOutliers (bufferB,&bufferBO,lastCountB,0.5);
	
	/* Show the results */
	avgA = avg (bufferAO,lastCountA);
	avgB = avg (bufferBO,lastCountB);
	minmaxA = getMinMax(bufferAO,lastCountA);
	minmaxB = getMinMax(bufferBO,lastCountB);
	

	printf ("\n\nA Count = %d",lastCountA);
	printf ("\nB Count = %d\n",lastCountB);
	printf ("\nA Min = %.6g", minmaxA.min);
	printf ("\nA Max = %.6g", minmaxA.max);
	printf ("\nB Min = %.6g", minmaxB.min);
	printf ("\nB Max = %.6g\n", minmaxB.max);
	printf ("\nAVG A = %.6g ",avgA);
	printf ("\nAVG B = %.6g \n",avgB);
	printf ("\nSD A = %.6g ",SDSamp(bufferA,lastCountA));
	printf ("\nSD B = %.6g \n",SDSamp(bufferB,lastCountB));



	printf ("\n*** Welch t-test Unpaired ***");

	p = PValueUnpaired (bufferAO,lastCountA,bufferBO,lastCountB);

	printf ("\nP-Value Two Sided = %.6g ",p);

	if (avgA <= avgB)
	{
		printf ("\nP-Value One Sided A < B = %.6g ", 0.5*p);
		printf ("\nP-Value One Sided A > B = %.6g ", 1 - 0.5*p);
	}
	else
	{
		printf ("\nP-Value One Sided A < B = %.6g ", 1 - 0.5*p);
		printf ("\nP-Value One Sided A > B = %.6g ", 0.5*p);
	}

	/* Clean up after ourselves */
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
