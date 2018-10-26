#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CSVLib.h>
#include <math.h>

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

double SD(double *buffer, int count)
{
	double mean, standardDeviation = 0.0;
	int i;

	mean = avg(buffer,count);

	for(i=0; i<count; ++i) standardDeviation += pow(buffer[i] - mean, 2);

	return sqrt(standardDeviation/count);
}

double PfromT (const double WELCH_T_STATISTIC, const double DEGREES_OF_FREEDOM)
{ 	
	const double a = DEGREES_OF_FREEDOM/2;
	double value = DEGREES_OF_FREEDOM/(WELCH_T_STATISTIC*WELCH_T_STATISTIC+DEGREES_OF_FREEDOM);
	
	if ((isinf(value) != 0) || (isnan(value) != 0)) return 1.0;

	printf ("\nT = %f",WELCH_T_STATISTIC);
	printf ("\ndf = %f\n",DEGREES_OF_FREEDOM);

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

double PValueUnpaired (double *ARRAY1, const size_t ARRAY1_SIZE, double *ARRAY2, const size_t ARRAY2_SIZE) 
{
	if (ARRAY1_SIZE <= 1) return 1.0;
	if (ARRAY2_SIZE <= 1) return 1.0;
	
	double fmean1 = 0.0, fmean2 = 0.0;
	
	fmean1 = avg(ARRAY1,ARRAY1_SIZE);
	fmean2 = avg(ARRAY2,ARRAY2_SIZE);
	
	if (fmean1 == fmean2) return 1.0;

	double unbiased_sample_variance1 = 0.0, unbiased_sample_variance2 = 0.0;

	for (size_t x = 0; x < ARRAY1_SIZE; x++) 
	{
		unbiased_sample_variance1 += (ARRAY1[x]-fmean1)*(ARRAY1[x]-fmean1);
	}

	for (size_t x = 0; x < ARRAY2_SIZE; x++) 
	{
		unbiased_sample_variance2 += (ARRAY2[x]-fmean2)*(ARRAY2[x]-fmean2);
	}
	
	unbiased_sample_variance1 = unbiased_sample_variance1/(ARRAY1_SIZE-1);
	unbiased_sample_variance2 = unbiased_sample_variance2/(ARRAY2_SIZE-1);

	const double WELCH_T_STATISTIC = (fmean1-fmean2)/sqrt(unbiased_sample_variance1/ARRAY1_SIZE+unbiased_sample_variance2/ARRAY2_SIZE);
	const double DEGREES_OF_FREEDOM = pow((unbiased_sample_variance1/ARRAY1_SIZE+unbiased_sample_variance2/ARRAY2_SIZE),2.0)/((unbiased_sample_variance1*unbiased_sample_variance1)/(ARRAY1_SIZE*ARRAY1_SIZE*(ARRAY1_SIZE-1))+(unbiased_sample_variance2*unbiased_sample_variance2)/(ARRAY2_SIZE*ARRAY2_SIZE*(ARRAY2_SIZE-1)));

	return PfromT(WELCH_T_STATISTIC,DEGREES_OF_FREEDOM);
}

double PValuePaired (double *ARRAY1, double *ARRAY2, const size_t ARRAY_SIZE) 
{
	double *ABdiff = NULL;
	double mean;
	double std;
	double WELCH_T_STATISTIC;
	double DEGREES_OF_FREEDOM = ARRAY_SIZE - 1;

	int i;
	
	ABdiff = malloc (ARRAY_SIZE * sizeof(double)); 
	if (ABdiff == NULL) exit (EXIT_FAILURE);

	for (i=0;i<ARRAY_SIZE;i++) ABdiff[i] = ARRAY1[i] - ARRAY2[i];

	mean = avg(ABdiff,ARRAY_SIZE);

	std = SD(ABdiff,ARRAY_SIZE);

	WELCH_T_STATISTIC = mean/(std/sqrt(ARRAY_SIZE-1));

	free (ABdiff);

	return PfromT(WELCH_T_STATISTIC,DEGREES_OF_FREEDOM);
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
	

	printf ("\nA Count = %d",lastCountA);
	printf ("\nB Count = %d\n",lastCountB);
	printf ("\nA Min = %f", minmaxA.min);
	printf ("\nA Max = %f", minmaxA.max);
	printf ("\nB Min = %f", minmaxB.min);
	printf ("\nB Max = %f\n", minmaxB.max);
	printf ("\nAVG A = %f ",avgA);
	printf ("\nAVG B = %f \n",avgB);
	printf ("\nSD A = %f ",SD(bufferA,lastCountA));
	printf ("\nSD B = %f \n",SD(bufferB,lastCountB));


	printf ("\n*** Unpaired ***");

	p = PValueUnpaired (bufferA,lastCountA,bufferB,lastCountB);

	printf ("\nP-Value Two Sided = %f ",p);

	if (avgA <= avgB)
	{
		printf ("\nP-Value One Sided A < B = %f ", 0.5*p);
		printf ("\nP-Value One Sided A > B = %f ", 1 - 0.5*p);
	}
	else
	{
		printf ("\nP-Value One Sided A < B = %f ", 1 - 0.5*p);
		printf ("\nP-Value One Sided A > B = %f ", 0.5*p);
	}

	if (lastCountA == lastCountB)
	{
		printf ("\n\n*** Paired ***");
		
		p = PValuePaired(bufferA,bufferB,lastCountA);
	
		printf ("\nP-Value Two Sided = %f ",p);

		if (avgA <= avgB)
		{
			printf ("\nP-Value One Sided A < B = %f ", 0.5*p);
			printf ("\nP-Value One Sided A > B = %f ", 1 - 0.5*p);
		}
		else
		{
			printf ("\nP-Value One Sided A < B = %f ", 1 - 0.5*p);
			printf ("\nP-Value One Sided A > B = %f ", 0.5*p);
		}
	}
		
	/* Clean up after ourselves */
	free (bufferA);
	bufferA = NULL;

	free (bufferB);
	bufferB =  NULL;

	printf ("\n\n");

    /* And we are out of here! */
    return 0;
}
