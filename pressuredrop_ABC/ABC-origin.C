/* ABC algorithm coded using C programming language */

/* Artificial Bee Colony (ABC) is one of the most recently defined algorithms by Dervis Karaboga in 2005, 
motivated by the intelligent behavior of honey bees. */

/* Referance Papers*/

/*D. Karaboga, AN IDEA BASED ON HONEY BEE SWARM FOR NUMERICAL OPTIMIZATION,TECHNICAL REPORT-TR06, Erciyes University, Engineering Faculty, Computer Engineering Department 2005.*/

/*D. Karaboga, B. Basturk, A powerful and Efficient Algorithm for Numerical Function Optimization: Artificial Bee Colony (ABC) Algorithm, Journal of Global Optimization, Volume:39, Issue:3,pp:459-171, November 2007,ISSN:0925-5001 , doi: 10.1007/s10898-007-9149-x */

/*D. Karaboga, B. Basturk, On The Performance Of Artificial Bee Colony (ABC) Algorithm, Applied Soft Computing,Volume 8, Issue 1, January 2008, Pages 687-697. */

/*D. Karaboga, B. Akay, A Comparative Study of Artificial Bee Colony Algorithm,  Applied Mathematics and Computation, 214, 108-132, 2009. */

/*Copyright © 2009 Erciyes University, Intelligent Systems Research Group, The Dept. of Computer Engineering*/

/*Contact:
Dervis Karaboga (karaboga@erciyes.edu.tr )
Bahriye Basturk Akay (bahriye@erciyes.edu.tr)
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ncurses.h>
#include <time.h>
#include <string.h>


/* Control Parameters of ABC algorithm*/
#define NP 10 /* The number of colony size (employed bees+onlooker bees)*/
#define FoodNumber NP/2 /*The number of food sources equals the half of the colony size*/
#define limit 50  /*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
#define maxCycle 500 /*The number of cycles for foraging {a stopping criteria}*/

/* Problem specific variables*/
#define D 5 /*The number of parameters of the problem to be optimized*/
#define lb 13 /*lower bound of the parameters. */
#define ub 15 /*upper bound of the parameters. lb and ub can be defined as arrays for the problems of which parameters have different bounds*/


#define runtime 1  /*Algorithm can be run many times in order to see its robustness*/


double Foods[FoodNumber][D]; /*Foods is the population of food sources. Each row of Foods matrix is a vector holding D parameters to be optimized. The number of rows of Foods matrix equals to the FoodNumber*/
double f[FoodNumber];  /*f is a vector holding objective function values associated with food sources */
double sen[FoodNumber][D]={0};
double fitness[FoodNumber]; /*fitness is a vector holding fitness (quality) values associated with food sources*/
double trial[FoodNumber]; /*trial is a vector holding trial numbers through which solutions can not be improved*/
double prob[FoodNumber]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
double solution [D]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
double ObjValSol; /*Objective function value of new solution*/
double Objsen[5]={0};
double FitnessSol; /*Fitness value of new solution*/
int neighbour, param2change; /*param2change corrresponds to j, neighbour corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
double GlobalMin; /*Optimum solution obtained by ABC algorithm*/
double GlobalParams[D]; /*Parameters of the optimum solution*/
double GlobalMins[runtime]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/
double r; /*a random number in the range [0,1)*/

/*a function pointer returning double and taking a D-dimensional array as argument */
/*If your function takes additional arguments then change function pointer definition and lines calling "...=function(solution);" in the code*/
typedef double (*FunctionCallback)(double sol[D], double sensitivity[D]); 

/*benchmark functions */
double sphere(double sol[D]);
double Rosenbrock(double sol[D]);
double Griewank(double sol[D]);
double Rastrigin(double sol[D]);
double pressuredrop(double sol[D], double sensitivity[D]);

/*Write your own objective function name instead of sphere*/
FunctionCallback function = &pressuredrop;

/*Fitness function*/
double CalculateFitness(double fun)
 {
	 double result=0;
	 if(fun>=0)
	 {
		 result=1/(fun+1);
	 }
	 else
	 {
		 result=1+fabs(fun);
	 }
	 return result;
 }

/*The best food source is memorized*/
void MemorizeBestSource()
{
   int i,j;
    
	for(i=0;i<FoodNumber;i++)
	{
	if (f[i]<GlobalMin)
		{
        GlobalMin=f[i];
        for(j=0;j<D;j++)
           GlobalParams[j]=Foods[i][j];
        }
	}
 }

/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */
/* Counters of food sources are also initialized in this function*/
void init(int index)
{
   int j;
   for (j=0;j<D;j++)
		{
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        Foods[index][j]=r*(ub-lb)+lb;
		solution[j]=Foods[index][j];
		}
	f[index]=function(solution,sen[index]);
	fitness[index]=CalculateFitness(f[index]);
	trial[index]=0;
}

/*All food sources are initialized */
void initial()
{
	int i;
	for(i=0;i<FoodNumber;i++)
	{
	init(i);
	}
	GlobalMin=f[0];
    for(i=0;i<D;i++)
    GlobalParams[i]=Foods[0][i];


}

void SendEmployedBees()
{
  int i,j;
  /*Employed Bee Phase*/
   for (i=0;i<FoodNumber;i++)
        {
        /*The parameter to be changed is determined randomly*/
        r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        param2change=(int)(r*D);
        
        /*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        neighbour=(int)(r*FoodNumber);

        /*Randomly selected solution must be different from the solution i*/        
        while(neighbour==i)
        {
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        neighbour=(int)(r*FoodNumber);
        }
        for(j=0;j<D;j++)
        solution[j]=Foods[i][j];

        /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        solution[param2change]=Foods[i][param2change]+(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2;

        /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        if (solution[param2change]<lb)
           solution[param2change]=lb;
        if (solution[param2change]>ub)
           solution[param2change]=ub;
        ObjValSol=function(solution,Objsen);
        FitnessSol=CalculateFitness(ObjValSol);
        
        /*a greedy selection is applied between the current solution i and its mutant*/
        if (FitnessSol>fitness[i])
        {
        /*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
        trial[i]=0;
        for(j=0;j<D;j++)
        Foods[i][j]=solution[j];
        f[i]=ObjValSol;
	memcpy(sen[i],Objsen,sizeof(Objsen));
        fitness[i]=FitnessSol;
        }
        else
        {   /*if the solution i can not be improved, increase its trial counter*/
            trial[i]=trial[i]+1;
        }


        }

        /*end of employed bee phase*/

}

/* A food source is chosen with the probability which is proportioal to its quality*/
/*Different schemes can be used to calculate the probability values*/
/*For example prob(i)=fitness(i)/sum(fitness)*/
/*or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/
void CalculateProbabilities()
{
     int i;
     double maxfit;
     maxfit=fitness[0];
  for (i=1;i<FoodNumber;i++)
        {
           if (fitness[i]>maxfit)
           maxfit=fitness[i];
        }

 for (i=0;i<FoodNumber;i++)
        {
         prob[i]=(0.9*(fitness[i]/maxfit))+0.1;
        }

}

void SendOnlookerBees()
{

  int i,j,t;
  i=0;
  t=0;
  /*onlooker Bee Phase*/
  while(t<FoodNumber)
        {

        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        if(r<prob[i]) /*choose a food source depending on its probability to be chosen*/
        {        
        t++;
        
        /*The parameter to be changed is determined randomly*/
        r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        param2change=(int)(r*D);
        
        /*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        neighbour=(int)(r*FoodNumber);

        /*Randomly selected solution must be different from the solution i*/        
        while(neighbour==i)
        {
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        neighbour=(int)(r*FoodNumber);
        }
        for(j=0;j<D;j++)
        solution[j]=Foods[i][j];

        /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        solution[param2change]=Foods[i][param2change]-(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2;

        /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        if (solution[param2change]<lb)
           solution[param2change]=lb;
        if (solution[param2change]>ub)
           solution[param2change]=ub;
        ObjValSol=function(solution,Objsen);
        FitnessSol=CalculateFitness(ObjValSol);
        
        /*a greedy selection is applied between the current solution i and its mutant*/
        if (FitnessSol>fitness[i])
        {
        /*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
        trial[i]=0;
        for(j=0;j<D;j++)
        Foods[i][j]=solution[j];
        f[i]=ObjValSol;
	memcpy(sen[i],Objsen,sizeof(Objsen));
        fitness[i]=FitnessSol;
        }
        else
        {   /*if the solution i can not be improved, increase its trial counter*/
            trial[i]=trial[i]+1;
        }
        } /*if */
        i++;
        if (i==FoodNumber)
        i=0;
        }/*while*/

        /*end of onlooker bee phase     */
}

/*determine the food sources whose trial counter exceeds the "limit" value. In Basic ABC, only one scout is allowed to occur in each cycle*/
void SendScoutBees()
{
int maxtrialindex,i;
maxtrialindex=0;
for (i=1;i<FoodNumber;i++)
        {
         if (trial[i]>trial[maxtrialindex])
         maxtrialindex=i;
        }
if(trial[maxtrialindex]>=limit)
{
	init(maxtrialindex);
}
}


/*Main program of the ABC algorithm*/
int main()
{
int iter,run,j;
double mean;
mean=0;
srand(time(NULL));

for(run=0;run<runtime;run++)
{

initial();
MemorizeBestSource();
for (iter=0;iter<maxCycle;iter++)
    {
    SendEmployedBees();
    CalculateProbabilities();
    SendOnlookerBees();
    MemorizeBestSource();
    SendScoutBees();

FILE *fp;
fp=fopen("GlobalMin","a");
if (!fp)
return -1;
fprintf(fp,"%d. iter: %e %e %e %e %e %e \n",iter+1,GlobalMin,GlobalParams[0],GlobalParams[1],GlobalParams[2],GlobalParams[3],GlobalParams[4]);
fclose(fp);

    }
for(j=0;j<D;j++)
		{
			printf("GlobalParam[%d]: %e\n",j+1,GlobalParams[j]);
		}
printf("%d. run: %e \n",run+1,GlobalMin);

GlobalMins[run]=GlobalMin;
mean=mean+GlobalMin;
}
mean=mean/runtime;
printf("Means of %d runs: %e\n",runtime,mean);
getch();
}


double sphere(double sol[D])
{
int j;
double top=0;
for(j=0;j<D;j++)
{
top=top+sol[j]*sol[j];
}
return top;
}

double Rosenbrock(double sol[D])
{
int j;
double top=0;
for(j=0;j<D-1;j++)
{
top=top+100*pow((sol[j+1]-pow((sol[j]),(double)2)),(double)2)+pow((sol[j]-1),(double)2);
}
return top;
}

 double Griewank(double sol[D])
 {
	 int j;
	 double top1,top2,top;
	 top=0;
	 top1=0;
	 top2=1;
	 for(j=0;j<D;j++)
	 {
		 top1=top1+pow((sol[j]),(double)2);
		 top2=top2*cos((((sol[j])/sqrt((double)(j+1)))*M_PI)/180);

	 }	
	 top=(1/(double)4000)*top1-top2+1;
	 return top;
 }
double pressuredrop(double sol[D], double sensitivity[D])
{
int j;
double top=0;
char s[120][120]={0};
char t[120][120]={0};
int  tline=0;
int   stime;
FILE *fp;
fp=fopen("constant/polyMesh/blockMeshDict","r");
if (!fp)
return -1;
while(tline<120)
	{
	fgets(s[tline],120,fp);
	tline++;
	}
fclose(fp);
fp=fopen("constant/polyMesh/blockMeshDict","w");
for(tline=0;tline<120;tline++)
	{
	if(tline==38)
	fprintf(fp,"spline 0 1 ((14 %f -0.5)(22 %f -0.5)(30 %f -0.5)(38 %f -0.5)(46 %f -0.5))\n",sol[0],sol[1],sol[2],sol[3],sol[4]);
	else if(tline==39)
	fprintf(fp,"spline 4 5 ((14 %f 0.5)(22 %f 0.5)(30 %f 0.5)(38 %f 0.5)(46 %f 0.5))\n",sol[0],sol[1],sol[2],sol[3],sol[4]);
	else
	fprintf(fp,"%s",s[tline]);
	}
fclose(fp);

tline=0;
fp=fopen("system/controlDict","r");
while(tline<120)
	{
	fgets(t[tline],120,fp);
	tline++;
	}
fclose(fp);
fp=fopen("system/controlDict","w");
for(tline=0;tline<120;tline++)
	{
	if(tline==92)
	fprintf(fp,"(0.014 %f 0)\n",double(0.001*sol[0]-0.000001));
	else if(tline==93)
	fprintf(fp,"(0.022 %f 0)\n",double(0.001*sol[1]-0.000001));
	else if(tline==94)
	fprintf(fp,"(0.030 %f 0)\n",double(0.001*sol[2]-0.000001));
	else if(tline==95)
	fprintf(fp,"(0.038 %f 0)\n",double(0.001*sol[3]-0.000001));
	else if(tline==96)
	fprintf(fp,"(0.046 %f 0)\n",double(0.001*sol[4]-0.000001));
	else
	fprintf(fp,"%s",t[tline]);
	}
fclose(fp);

system("blockMesh");
system("rm -r 1000/");
system("rm -r 2000/");
system("simpleFoam");
system("lighthillCostFunction > cost");
/*
system("cp 0/Ua  1000/");
system("cp 0/pa  1000/");
system("sed -i s/'\\(endTime[ \t]*\\) 1000;'/'\\1 2000;'/g system/controlDict");
system("sed -i s/'0.5'/'0.15'/g system/fvSolution");
system("sed -i s/'0.75'/'0.2'/g system/fvSolution");
system("adjointSensitivityFoam");
system("sed -i s/'\\(endTime[ \t]*\\) 2000;'/'\\1 1000;'/g system/controlDict");
system("sed -i s/'0.15'/'0.5'/g system/fvSolution");
system("sed -i s/'0.2'/'0.75'/g system/fvSolution");
*/
fp=fopen("cost","r");
for(tline=0;tline<120;tline++)
	{
	if(tline==42)
	fscanf(fp,"%lf",&top);
	else
	fgets(s[tline],120,fp);
	}
fclose(fp);

fp=fopen("postProcessing/probes/1000/sensitivity","r");
for(tline=0;tline<6;tline++)
	{
	if(tline==4)
	fscanf(fp,"%d%lf%lf%lf%lf%lf",&stime,&sensitivity[0],&sensitivity[1],&sensitivity[2],&sensitivity[3],&sensitivity[4]);
	else

	fgets(s[tline],120,fp);

	}
fclose(fp);
return top;
}


 double Rastrigin(double sol[D])
 {
	 int j;
	 double top=0;

	 for(j=0;j<D;j++)
	 {
		 top=top+(pow(sol[j],(double)2)-10*cos(2*M_PI*sol[j])+10);
	 }
	 return top;
 }
