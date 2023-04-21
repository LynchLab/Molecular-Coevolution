/* 

This program estimates the quasi-steady-state features of a two-locus / two-allele system in a finite population, allowing for arbitrary selection and reversible mutation.

The goal is to estimate the steady-state distribution of the four genotypes over time, and the evolutionary substitution rate per site relative to the neutral expectation. 

The equilibrium long-term result is a balance between the forces of mutation, selection, and random genetic drift. 

Results are obtained for an array of actual population sizes (N), which are constant in time.

The population experiences sequential episodes of mutation, recombination, selection, and random genetic drit.

Haploidy is assumed.

Haplotype designations: AB = [1][1], Ab = [1][0], aB = [0][1], and ab = [0][0], with fitness schemes 1, 1 - sb, 1 - sa, and 1 + sab.

All mutation, recombination, and selection processes are treated deterministically, prior to random sampling.

Mutation:
	There are just two types of mutations at each type of site: - to +, ux01 (beneficial); and + to -, ux10 (deleterious), where x = a or b; the rates are allowed to differ among loci.
	
Drift:
	Separate effective population sizes are allowed for the two loci, but if complete linkage is assumed, these should be set equal to each other (i.e., set neratio = 1.0). 

	The two loci are allowed to have different Ne, which requires application of two drift samplings, the first at Ne = Nea over the full genotype; and the second applying Nex only to the B locus, 
	keeping the previous A frequencies constant. 

	Letting 1 - (1/Neb) = [1 - (1/Nea)] * [1 - (1/Nex)], Nea >= Neb, leads to 

	Nex = Neb * (Nea - 1) / (Nea - Neb) to be used in the 2nd interval of sampling.
	


Running averages of the population features are kept track of: equilibrium frequencies of the four types; average Ne at the two loci.

The run starts with an even genotype-frequency distribution, but this can be modified internally.

After a burnin, statistics are then taken at intervals of N/xinc generations, with a total of ngen sampling intervals.

An array of population sizes to run is set internally.


*/



/* ********************************************************************************************************************** */

#define ua		0.000000001					/* mutation rate of A to a */

#define ub		0.0000001					/* mutation rate of B to b */

#define muta		1.0						/* ratio of beneficial to deleterious mutation rates, locus A */

#define mutb		1.0						/* ratio of beneficial to deleterious mutation rates, locus B */

#define scoa		0.000001				    /* selection coefficient against aB (positive if deleterious) */

#define scob		0.00001				    /* selection coefficient against Ab (positive if deleterious) */

#define scoab		0.00000					/* selection coefficient in favor of AB (positive if advantageous) */

#define neratio		0.01					/* ratio of Ne at locus B to that at locus A (array in program is for locus A, so this should be <= 1.0, with B having the smaller Ne, as in an organelle genome) */

#define recr		0.5						/* recombination rate: recr = 0.5 is free recombination. */

#define xinc		25						/* statistics to be recorded every ne/xinc generations */

#define burnin		50000					/* number of initial burn-in sampling increments, ignored in the statistics */

#define tintprint	1000000					/* number of sampling increments between screen printing */



#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>
#include    <string.h> 
#include    <gsl/gsl_rng.h>
#include    <gsl/gsl_randist.h>



/* ********************************************************************************************************** */


/* Point to the output file. */

FILE *stream;
char filename[100];

 void main(int argc, char *argv[]) 
{                                                                               
 int f0, f1;                                                                     
                                                                                 
 if (argc > 1)                                                                 
 {                                                                           
     f0 = atoi(argv[1]);                                                             
     f1=f0;                                                                      
     sprintf(filename, "dataouts65_%d.txt", f0);                                        
 }                                                                               
 else                                                                          
 {                                                                           
     f0 = 1;                                                                     
     f1 = 19;                                                                    
     sprintf(filename, "dataouts65.txt");                                               
 }  



/* ********************************************************************************************************** */


/* Call for the random number generator. */

 static gsl_rng* rand_new;
 static gsl_rng_type * T;
 gsl_rng_env_setup();
 if (!rand_new)
 {
	 rand_new = gsl_rng_alloc(gsl_rng_taus2);
	 gsl_rng_set(rand_new, time(NULL) + (20*f0));
 }



/* ***************************************************************************************** */

	/* MAIN BODY OF PROGRAM. */


	int iga, igb;										/* counters for the classes */

	long igen;											/* generation counter */

	double tint;										/* interval number for statistics */

	long nea, neb, nex;									/* effective population sizes */

	double u10a, u01a, u10b, u01b;						/* mutation rates; 01 = beneficial, 10 = deleterious */

	double selcoa, selcob, selcoab;					    /* selection coefficients */

    double crate, reccor;                               /* potentially scaled recombination rate */

	int kfacn, kfacs, kfac;							    /* scaling factor for speeding up runs, from scalef[] */
    double kfacsdbl, kfacrec;

	long efpopn[40];									/* effective population sizes at locus A to run -- NOTE THESE ARE SCALED DOWN ACCORDING TO SCALEF TO INCREASE RUN SPEED */
    double rlng[40];                                    /* SCALING FACTORS FOR SPEEDING UP RUNS */

    double ngens;                                       /* time iterations in run */

	int itera;											/* counter for population-size iterations */

	long increment;										/* increment between surveys */
	long tcount;										/* counter for printing to screen to monitor simulation progress */
	long counter;										/* initial number of surveys to be skipped */

	double meanfit;									    /* mean fitness */

	double ldiseq;										/* linkage disequilibrium */

	double wfit[2][2];									/* genotypic fitnesses */

	double p0[2][2];									/* genotypic frequencies at start of generation */

	double psel[2][2];									/* after selection */

	double pmutm[2][2];								    /* after mutation */

	double prec[2][2];									/* after recombination */

	double pgtypexp[2][2];								/* expected frequencies prior to first epsisode of random drift */

	double pnew[2][2];									/* frequencies after 2nd sampling episode */

	double pa1tot, pa0tot;								/* conditional frequencies of A/a alleles */

	double sumfreqab[2][2];							    /* summations of genotype frequencies*/

	double totw;										/* summations for grand means and variances */

	double sump;										/* sum of frequencies */

	double pp;                                         	/* probability associated with the binomial for drift */
	
	long ntot;											/* integer associated with the binomial for drift */
	long draw;											/* drift drawn from the binomial */
	
	
	double epoi, rnum;									/* terms for Poisson draws */

	double meanw, grandmeanw;							/* generational mean for fitness */

	double mean[3][3];									/* mean genotypic frequencies */

	int oldfix, newfix;									/* indicators for fixation states */

	double fixgens, numfix[5][5], genfix[5][5];		    /* counters for numbers and times of fixations of different types */
	
	double totgens, totratea, totrateb;

	double meanfixtime[5][5];							/* mean fixation times */
	double fixrate[5][5];								/* mean fixation rates */

	double ratea, rateb, neutratea, neutrateb;			/* overall average observed substitution rates and expected netural rates */

    double start, stop, time;                                                   


	/* Open the output file. */

	remove("dataouts65.txt ");





	/* Population sizes to run. */

	efpopn[19] = 1000000000;
	efpopn[18] = 465000000;
	efpopn[17] = 216000000;
	efpopn[16] = 100000000;
	efpopn[15] = 46500000;
	efpopn[14] = 21600000;
	efpopn[13] = 10000000;
	efpopn[12] = 4650000;
	efpopn[11] = 2160000;
	efpopn[10] = 1000000;
	efpopn[9] = 465000;
	efpopn[8] = 216000;
	efpopn[7] = 100000;
	efpopn[6] = 46500;
	efpopn[5] = 21600;
	efpopn[4] = 10000;
	efpopn[3] = 4650;
	efpopn[2] = 2160;
	efpopn[1] = 1000;




	/* Number of sampling increments in run; each increment is (ne/10) generations */

	rlng[19] = 50000000.0;
	rlng[18] = 50000000.0;
	rlng[17] = 100000000.0;
	rlng[16] = 100000000.0;
	rlng[15] = 400000000.0;
	rlng[14] = 400000000.0;
	rlng[13] = 1000000000.0;
	rlng[12] = 2000000000.0;
	rlng[11] = 2000000000.0;
	rlng[10] = 2000000000.0;
	rlng[9] = 2000000000.0;
	rlng[8] = 2000000000.0;
	rlng[7] = 2000000000.0;
	rlng[6] = 4000000000.0;
	rlng[5] = 10000000000.0;
	rlng[4] = 10000000000.0;
	rlng[3] = 20000000000.0;
	rlng[2] = 20000000000.0;
	rlng[1] = 20000000000.0;




for (itera = f0; itera <= f1; ++itera){							/* Start iterations over the set of population sizes and mutation rates. */
        
        stream=fopen(filename, "a");

        /* Set the run length. */
        
        ngens = rlng[itera];



		/* Set the initial population-genetic parameters. */

		nea = efpopn[itera];										/* effective population size for locus A */
		neb = int(((double) nea) * neratio);						/* effective population size for locus B */

		nex = 0;
		if (neb != nea) {
			nex = int(((double)neb) * (((double)nea) - 1.0) / (((double)nea) - ((double)neb))); }		/* effective population size to be used for second sampling episode to allow for lower Ne at locus B */
																										/* with the double sampling scheme used for locus B, this renders the desired effective size = neb */


    /* Scale the effective population size down, and u, s, and c up to speed up runs. */


        kfacn = neb / 1000;
    
        if ( (scoa*scob) == 0.0 ) {
            kfacsdbl = ((double) kfacn);}
        else { kfacsdbl = 0.1 / scoa;}                          /* make sure the scaled selection coefficient does not exceed 0.1 */ 
        
        kfacs = int(kfacsdbl);

        if (kfacn < 1.0) {
            kfacn = 1.0; }
        if (kfacs < 1.0) {
            kfacs = 1.0; }

        if ((kfacn < kfacs) ) {
            kfac = kfacn; }
        else {kfac = kfacs;}



        reccor = recr;
        
        crate = ((double) kfac) * recr;                         /* make sure the scaled recombination rate does not exceed 0.5 */
        
        if ( (recr > 0.0) && (recr < 0.5) ){
        
            if (crate >= 0.5) {
                kfacrec = 0.5 / recr;
                kfac = int(kfacrec);}
            
            if (recr < 0.5) {
                reccor = ((double) kfac) * recr; }
        }


		nea = nea / kfac;
		neb = neb / kfac;
		nex = nex / kfac;

		u10a = ((double) kfac) * ua;							/* deleterious and beneficial mutation rates at locus A */
		u01a = muta * u10a;

		u10b = ((double) kfac) * ub;							/* deleterious and beneficial mutation rates at locus B */
		u01b = mutb * u10b;

		selcoa = ((double) kfac) * scoa;						/* selection coefficients */
		selcob = ((double) kfac) * scob;
		selcoab = ((double) kfac) * scoab;




		wfit[1][1] = 1.0;										/* genotypic fitnesses */
		wfit[0][1] = 1.0 - selcoa;
		wfit[1][0] = 1.0 - selcob;
		wfit[0][0] = 1.0 + selcoab;

		p0[1][1] = 0.25;										/* set the initial genotype frequencies to be random. THIS IS ARBITRARY, AND THE HISTORICAL EFFECT ON RUNS IS MINIMAL WITH A LONG BURN-IN PERIOD. */
		p0[1][0] = 0.25;
		p0[0][1] = 0.25;
		p0[0][0] = 0.25;


		/* Initiate the counters for estimating transition rates between states. */

		oldfix = 1;												
		fixgens = 0.0;
		totgens = 0.0;
		totratea = 0.0;
		totrateb = 0.0;

		for (iga = 1; iga <= 4; ++iga) {
			for (igb = 1; igb <= 4; ++igb) {
				numfix[iga][igb] = 0.0;
				genfix[iga][igb] = 0.0;
				meanfixtime[iga][igb] = 0.0;
				fixrate[iga][igb] = 0.0; } }


		/* Initiate the genotype frequencies and counters. */

		for (iga = 0; iga <= 1; ++iga) {
			for (igb = 0; igb <= 1; ++igb) {
				pmutm[iga][igb] = 0.0;									/* zero the various allele-frequency counters */
				psel[iga][igb] = 0.0;
				pgtypexp[iga][igb] = 0.0;
				pnew[iga][igb] = 0.0;
				sumfreqab[iga][igb] = 0.0; 	} }

		igen = 0;
		tcount = 0;
		tint = 0.0;
		counter = 0;
		totw = 0.0;

		increment = nea / xinc;											/* increment in generations between statistic calculations (set as a fraction of nea). */




		/* ******************************************************************************************************************************************* */


		/* Iterate the recursion equations to obtain the equilibrium expectations. */

		while (tint < ngens)  										/* iterate until the stopping criterion has been met. */
		{
			igen = igen + 1;


			/* Impose selection on the genotypic classes. */

			meanfit = 0.0;

			for (iga = 0; iga <= 1; ++iga) {						/* calculate mean relative fitness */
				for (igb = 0; igb <= 1; ++igb) {
					meanfit = meanfit + (p0[iga][igb] * wfit[iga][igb]); } 	}

			for (iga = 0; iga <= 1; ++iga) {						/* genotypic frequencies after selection */
				for (igb = 0; igb <= 1; ++igb) {
					psel[iga][igb] = p0[iga][igb] * wfit[iga][igb] / meanfit; } }




			/* Impose mutation on the post-selection genotypic classes. */

			pmutm[1][1] = (u01a * u01b * psel[0][0])                 + (u01a * (1.0 - u10b) * psel[0][1])         + ((1.0 - u10a) * u01b * psel[1][0])         + ((1.0 - u10a) * (1.0 - u10b) * psel[1][1]);

			pmutm[1][0] = (u01a * (1.0 - u01b) * psel[0][0])         + (u01a * u10b * psel[0][1])                 + ((1.0 - u10a) * (1.0 - u01b) * psel[1][0]) + ((1.0 - u10a) * u10b * psel[1][1]);

			pmutm[0][1] = ((1.0 - u01a) * u01b * psel[0][0])         + ((1.0 - u01a) * (1.0 - u10b) * psel[0][1]) + (u10a * u01b * psel[1][0])                 + (u10a * (1.0 - u10b) * psel[1][1]);

			pmutm[0][0] = ((1.0 - u01a) * (1.0 - u01b) * psel[0][0]) + ((1.0 - u01a) * u10b * psel[0][1])         + (u10a * (1.0 - u01b) * psel[1][0])         + (u10a * u10b * psel[1][1]);




			/* Impose recombination on the post-mutation classes. */

			ldiseq = (pmutm[0][1] * pmutm[1][0]) - (pmutm[1][1] * pmutm[0][0]);

			prec[1][1] = pmutm[1][1] + (reccor * ldiseq);
			prec[0][0] = pmutm[0][0] + (reccor * ldiseq);
			prec[0][1] = pmutm[0][1] - (reccor * ldiseq);
			prec[1][0] = pmutm[1][0] - (reccor * ldiseq);





			/* Reset the next generation's expected genotype frequencies, and ensure that they sum to 1.0. */

			sump = prec[1][1] + prec[0][1] + prec[1][0] + prec[0][0];

			for (iga = 0; iga <= 1; ++iga) {
				for (igb = 0; igb <= 1; ++igb) {
					pgtypexp[iga][igb] = prec[iga][igb] / sump; 
					p0[iga][igb] = 0.0;	}}
				
				


			/* Sample the population for new genotype frequencies after the first two-locus sampling episode. */
					/* This allows for sampling at the first level with nea, and is all that is required if nea = neb. */
			
			ntot = nea;
			sump = 0.0;

			for (iga = 0; iga <= 1; ++iga) {
				for (igb = 0; igb <= 1; ++igb) {

					if ((pgtypexp[iga][igb] > 0.0) && (ntot > 0) && (sump < 1.0))  {
						pp = pgtypexp[iga][igb] / (1.0 - sump);												/* this is the remaining frequency to sample */

						if (pp >= 1.0000000000000) {														
							draw = ntot;
							p0[iga][igb] = ((double)draw) / ((double)nea);
						}

						else { 
							draw = gsl_ran_binomial_tpe(rand_new, pp, ntot);
							p0[iga][igb] = ((double)draw) / ((double)nea);
						}

						ntot = ntot - draw;
						sump = sump + pgtypexp[iga][igb];

					}
				}
			}



			/* Sample the population for new genotype frequencies after the second B-locus sampling episode, if neb <> nea, using nex as the effective population size. */

			if (neratio != 1.0) {

				for (iga = 0; iga <= 1; ++iga) {
					for (igb = 0; igb <= 1; ++igb) {
						pnew[iga][igb] = 0.0; } }

				pa1tot = p0[1][1] + p0[1][0];									/* frequency of allele A at first locus */
				pa0tot = p0[0][1] + p0[0][0];									/* frequency of allele a at first locus */

				if (pa1tot > 0.0) {

					pp = p0[1][1] / pa1tot;
					ntot = int(pa1tot * ((double)nex));
					if (ntot < 1){
						ntot = 1; }

					if (pp == 1.0) {
						pnew[1][1] = p0[1][1]; }
					else if (pp == 0.0) {
						pnew[1][0] = p0[1][0]; 	}
					else {
						draw = gsl_ran_binomial_tpe(rand_new, pp, ntot);
						pnew[1][1] = pa1tot * ((double)draw) / ((double)ntot);
						pnew[1][0] = pa1tot - pnew[1][1];	}
				}

				if (pa0tot > 0.0) {

					pp = p0[0][0] / pa0tot;
					ntot = int(pa0tot * ((double)nex));
					if (ntot < 1){
						ntot = 1; }

					if (pp == 1.0) {
						pnew[0][0] = p0[0][0]; }
					else if (pp == 0.0) {
						pnew[0][1] = p0[0][1]; 	}
					else {
						draw = gsl_ran_binomial_tpe(rand_new, pp, ntot);
						pnew[0][0] = pa0tot * ((double)draw) / ((double)ntot);
						pnew[0][1] = pa0tot - pnew[0][0]; 	}
				}

				for (iga = 0; iga <= 1; ++iga) {
					for (igb = 0; igb <= 1; ++igb) {
						p0[iga][igb] = pnew[iga][igb];
						if (p0[iga][igb] < 0.0){
							p0[iga][igb] = 0.0;	}}}
			}




			/* Check for new fixation, and if it occurs, record statistics; cutoff frequency is arbitrarily set as 0.999. */
					/* Legend: AB = 1; Ab = 2; aB = 3; ab = 4. */
					/* Note that when N*u > 0.01, fixations in the above sense often do not occur, owing to mutation pressure. */ 

			totgens = totgens + 1.0;

			if (p0[1][1] > 0.999){ newfix = 1; }
			else if (p0[1][0] > 0.999){ newfix = 2; }
			else if (p0[0][1] > 0.999){ newfix = 3; }
			else if (p0[0][0] > 0.999){ newfix = 4; }
			else { newfix = 0; }

			fixgens = fixgens + 1.0;

			if ((newfix != oldfix) && (newfix != 0)) {
				numfix[oldfix][newfix] = numfix[oldfix][newfix] + 1.0;
				genfix[oldfix][newfix] = genfix[oldfix][newfix] + fixgens;

				if ((oldfix==1 && newfix==3) || (oldfix==1 && newfix==4) || (oldfix==2 && newfix==3) || (oldfix==2 && newfix==4) || (oldfix==3 && newfix==1) || (oldfix==3 && newfix==2) || (oldfix==4 && newfix==1) || (oldfix==4 && newfix==2)) {
					totratea = totratea + 1.0; }
					
				if ((oldfix==1 && newfix==2) || (oldfix==1 && newfix==4) || (oldfix==2 && newfix==3) || (oldfix==2 && newfix==1) || (oldfix==3 && newfix==4) || (oldfix==3 && newfix==2) || (oldfix==4 && newfix==1) || (oldfix==4 && newfix==3)) {
					totrateb = totrateb + 1.0; }

				oldfix = newfix;
				fixgens = 0.0; 	
			}



			/* Calculate the summary statistics if the sampling interval is completed. */

			if (igen == increment) {
				igen = 0;
				counter = counter + 1;

				if (counter > burnin) {

					meanw = 0.0;

					for (iga = 0; iga <= 1; ++iga) {
						for (igb = 0; igb <= 1; ++igb) {
							meanw = meanw + (p0[iga][igb] * wfit[iga][igb]);

							sumfreqab[iga][igb] = sumfreqab[iga][igb] + p0[iga][igb];} 	}

					totw = totw + meanw;
					tint = tint + 1.0;
					tcount = tcount + 1;


					if (tcount > tintprint) {

						for (iga = 0; iga <= 1; ++iga) {									/* cumulative mean genotypic frequencies */
							for (igb = 0; igb <= 1; ++igb) {
								mean[iga][igb] = sumfreqab[iga][igb] / tint; } }

						printf("%9d, %9d, %9d, %9d, %10.0f, %9.5f, %9.5f, %6.5f, %6.5f, %6.5f, %9.5f, %6.5f \n", (nea*kfac), (neb*kfac), (nex*kfac), kfac, tint, (totw / tint),
							mean[1][1], mean[1][0], mean[0][1], mean[0][0], (mean[1][1] + mean[1][0]), (mean[1][1] + mean[0][1]));

						tcount = 0;
					}
				}
			}								/* ends the summary statistic analysis for this point */
		}									/* ends the loop for generations for this population size. */



		/* Calculate the final statistics. */
		
		grandmeanw = totw / tint;								/* mean fitness */

		for (iga = 1; iga <= 4; ++iga) {										/* mean transtion times from one fixed state to another */
			for (igb = 1; igb <= 4; ++igb) {
				if (numfix[iga][igb] >= 1.0) {
					meanfixtime[iga][igb] = genfix[iga][igb] / numfix[iga][igb]; } }	}

		for (iga = 1; iga <= 4; ++iga) {										/* transition rates */
			for (igb = 1; igb <= 4; ++igb) {
				if (iga != igb) {
					if (meanfixtime[iga][igb] > 0.0) {
						fixrate[iga][igb] = 1.0 / meanfixtime[iga][igb];	} }	} }


		/* Rates of substitution at the two loci. */

		ratea = totratea / totgens;
		rateb = totrateb / totgens;

		neutratea = 2.0 * u10a * u01a / (u10a + u01a);							/* neutral rates of substitution */
		neutrateb = 2.0 * u10b * u01b / (u10b + u01b);



		/* Print the output. */

		fprintf(stream, " %11d, %11d, %11d ,, %12.11f, %12.11f, %4.3f, %4.3f ,, %12.11f,, %12.11f, %12.11f, %12.11f ,,  %6d ,, %17.0f, %13d, %17.0f ,, %12.11f ,, %12.11f, %12.11f ,, %12.11f, %12.11f, %12.11f, %12.11f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f,, %12.11f, %9.1f ,, %12.11f, %12.11f ,, %12.11f, %12.11f ,, %12.11f, %12.11f \n  ",
			(nea*kfac), (neb*kfac), (nex*kfac),
			ua, ub, muta, mutb, recr, scoa, scob, scoab,
			kfac, ngens, burnin, (tint*((double)increment)),
			grandmeanw,
			(mean[1][1] + mean[1][0]), (mean[1][1] + mean[0][1]), 
			mean[1][1], mean[1][0], mean[0][1], mean[0][0],
			fixrate[1][2], numfix[1][2],
			fixrate[1][3], numfix[1][3],
			fixrate[1][4], numfix[1][4],
			fixrate[2][1], numfix[2][1],
			fixrate[2][3], numfix[2][3],
			fixrate[2][4], numfix[2][4],
			fixrate[3][1], numfix[3][1],
			fixrate[3][2], numfix[3][2],
			fixrate[3][4], numfix[3][4],
			fixrate[4][1], numfix[4][1],
			fixrate[4][2], numfix[4][2],
			fixrate[4][3], numfix[4][3],
			ratea, rateb,
			neutratea, neutrateb,
			(ratea / neutratea), (rateb / neutrateb));


		printf("\n");

		fclose(stream);


	}									/* End of the population-size loop. */


exit(0);

}





