//Pso Code For Many Objective Problem
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include<math.h>
#define URAND	((double)rand()/((double)RAND_MAX))

#define INITRAND srand(time(0))
 
//all global variables 
	int k, jrand, numofFE = 0;
	
	int r1, r2, r3;
  	double 	Xl= 0,Xu=1;
   	int NPPSO,NOPSO, GmaxPSO, c, *indexPSO , s = 1;
	double c1=1.0,c2=1.0,W=0.4; 
	int ***populPSO, ***nextPSO, ***ptrPSO, *iptrPSO,***velocity, *UPSO,**pbest,*gbest=5000000000,*arr,*arrUI,*arrb;
	int NP, Gmax, NO,c, *index,DPSO,D;
   	int  ***popul, ***next, ***ptr, *iptr, *U;
	double CR = 0.9, F = 0.9,min_value = DBL_MAX, totaltime = 0.0;
  	   clock_t starttime, endtime;
		int countPSO=0;
int count=0;
//Vrp data 
int NumberOfVechiles=0;
int routes[25][50];
int citysToRoute[25];
int VrootCapacity[25];
double weightCector[5]={0.5,0.2,0.1,0.1,0.1};
int bc[1000],ec[1000];
int city1,city2;
 long int** mat;
int totalcost=0,j;
  int n=101;
 int x1,x2,x,y2,y;
  int yy ;
double AA;
  unsigned int xy;
//double dist;
FILE *fid,*ipdata;
FILE *fod;
double sum=0;



 int main()
{


	time_t mytime;
		mytime = time(NULL);
	//char *aa=strcat("TEST/",ctime(&mytime));
	char *aa=ctime(&mytime);
	printf("%s",aa);
	fod = fopen(aa,"w");
	int i,j,k;
	//taking iNPPSOut from user
   	/*printf("Entere the population size::\n");
	scanf("%d",&NPPSO);
	printf("Entere the Number of Dimensions::\n");
	scanf("%d",&DPSO);
	printf("Entere the Number of Iterations::\n");
	scanf("%d",&GmaxPSO);
	printf("Entere the Number of Objectives::\n");
	scanf("%d",&NOPSO);*/
	NPPSO=10;
	DPSO=100;
	GmaxPSO=200;
	NOPSO=5;
	NP=NPPSO;
	D=DPSO;
	NO=NOPSO;
	Gmax=5000;
	fprintf(fod,"\nMulti objective HDEPSO Algorithm\n");
	fprintf(fod,"\nPopulation Size::::%d",NP);
	fprintf(fod,"\nNumber of Dimensions::::%d",D);
	fprintf(fod,"\nNumber of Iterations::::%d",Gmax);
	//Memory Allocation to Pointers
	populPSO = (int ***)malloc(NPPSO*NOPSO*sizeof(int **));
  		
	velocity = (int ***)malloc(NPPSO*NOPSO*sizeof(int *));
  		 if (velocity == NULL) perror("malloc");
	nextPSO = (int ***)malloc(NPPSO*NOPSO*sizeof(int **));
   		if (nextPSO == NULL) perror("malloc");
	 UPSO = (int *)malloc((DPSO+1)*sizeof(int));
  		 if (UPSO == NULL) perror("malloc");
	 gbest = (int *)malloc(NPPSO*sizeof(int));
  		 if (gbest== NULL) perror("malloc");
	 iptrPSO = (int *)malloc((DPSO+1)*sizeof(int));
  		 if (iptrPSO== NULL) perror("malloc");
	 indexPSO = (int *)malloc(NOPSO*sizeof(int));
  		 if (indexPSO == NULL) perror("malloc");
	pbest = (int **)malloc(NPPSO*sizeof( int));
  		 if (pbest == NULL) perror("malloc");

	popul = (int ***)malloc(NP*NO*sizeof(int **));
  		 if (popul == NULL) perror("malloc");
	next = (int ***)malloc(NP*NO*sizeof(int **));
   		if (next == NULL) perror("malloc");
	 U = (int *)malloc((D+1)*sizeof(int));
  		 if (U == NULL) perror("malloc");
	 iptr = (int *)malloc((D+1)*sizeof(int));
  		 if (U == NULL) perror("malloc");
	 index = (double *)malloc(NO*sizeof(double));
  		 if (index == NULL) perror("malloc");
		arr = (int *)malloc((D+1)*sizeof(int));



arr = (int *)malloc((DPSO+1)*sizeof(int));
  	 if (arr == NULL) perror("malloc");
//arrb = (int *)malloc((D+1)*sizeof(int));
  //	 if (arrb == NULL) perror("malloc");
	arrUI= (int *)malloc((DPSO+1)*sizeof(int));
  	 if (arrUI == NULL) perror("malloc");
//printf("Array");
//generating array of 0 to D
for(i=0;i<DPSO;i++)
	{
	arr[i]=i+1;
//fscanf(ipdata,"%d",&arrb[i]);
	//arr[i]=arrb[i];
        //printf("%d\t",arr[i]);
	//fprintf(fid,"%d\t",arr[i]);
	}
//shffle array by swapping
void swap(int *a,int *b)
{
int temp=*a;
*a=*b;
*b=temp;
}
//get the random indexes
void randm(int arr[],int n)
{
int i;
srand(time(NULL));
for(i=0;i<n;i++)
	{
int j=rand()%(i+1);
swap(&arr[i],&arr[j]);
	}
}

FILE *file;
		mat=malloc(10000*sizeof(int));
		for(i=0;i<n;++i)
		mat[i]=malloc(8*sizeof(int));

 		file=fopen("cdp101.txt", "r");
		if (file==NULL)
		{
			printf("\n Error: Cannot open input file for reading \n");
		}
		 for(i = 0; i < n; i++)
 		 {
  		   for(j = 0; j < 8; j++)
      			{
      			 if (!fscanf(file, " %d", &mat[i][j]))
			   break;
     // mat[i][j] -= '0'; 
      printf("%d \t",mat[i][j]);
			bc[i]=mat[i][5];
      			  ec[i]=mat[i][6]+30;
      			}
	 printf("\n");
  		}



	starttime = clock();
	printf("Initial Population is::\n");
		for(i=0;i<DPSO;i++)
		{
		printf("Dim%d\t\t",(i+1));
		}
		printf("\n");
	   	for (k=0; k < NOPSO; k++)
		{
			printf("\NOPSObjective Function NOPSO::%d\n",(k+1));
			populPSO[k] = (int *)malloc(NPPSO*(DPSO+1)*sizeof(int));
			    //  if (popul[k][i] == NULL) perror("malloc");
			velocity[k] =(int *)malloc(NPPSO*(DPSO+1)*sizeof(int));
		        if (velocity[k] == NULL) perror("malloc");
			nextPSO[k] = (int *)malloc(NPPSO*(DPSO+1)*sizeof(int));
			      if (nextPSO[k] == NULL) perror("malloc");
			pbest[k] = (int *)malloc(NPPSO*(DPSO+1)*sizeof(int));
			      if (pbest[k] == NULL) perror("malloc");
			for (i=0; i < NPPSO; i++)
			{
			      populPSO[k][i] = (int *)malloc((DPSO+1)*sizeof(int));
			      if (populPSO[k][i] == NULL) perror("malloc");
				velocity[k][i] = (int *)malloc((DPSO+1)*sizeof(int));
			      if (velocity[k][i] == NULL) perror("malloc");
				randm(arr,DPSO);
			      for (j=0; j < DPSO; j++)
				{
					populPSO[k][i][j] = arr[j];
					printf("%d\t",populPSO[k][i][j]);
					velocity[k][i][j]=arr[j];;
				}
				if(k==0)
				{
			     	populPSO[k][i][DPSO] = Distance(populPSO[k][i],DPSO);
				printf("%d\t",populPSO[k][i][DPSO]);
				printf("\n");
				}
				if(k==1)
				{
			     	populPSO[k][i][DPSO] = NumberOfV(populPSO[k][i],DPSO);
				//popul[k][i][D] = func(popul[k][i]);
				printf("%d\t",populPSO[k][i][DPSO]);
				printf("\n");
				}
				if(k==2)
				{
			     	populPSO[k][i][DPSO] = Makespan(populPSO[k][i],DPSO);
				//popul[k][i][D] = func(popul[k][i]);
				printf("%d\t",populPSO[k][i][DPSO]);
				printf("\n");
				}
				if(k==3)
				{
			     	populPSO[k][i][DPSO] = WaitingTime(populPSO[k][i],DPSO);
				//popul[k][i][D] = func(popul[k][i]);
				printf("%d\t",populPSO[k][i][DPSO]);
				printf("\n");
				}
				if(k==4)
				{
			     	populPSO[k][i][DPSO] = TotalDelay(populPSO[k][i],DPSO);
				//popul[k][i][D] = func(popul[k][i]);
				printf("%d\t",populPSO[k][i][DPSO]);
				printf("\n");
				}
			//printf("\n");
			     	numofFE++;

			     nextPSO[k][i] = (int *)malloc((DPSO+1)*sizeof(int));
			      if (nextPSO[k][i] == NULL) perror("malloc");
			}//i
		}//k
for (j=0;j<GmaxPSO;j++)
	{
		
		findPbest(j);
		findGbest();
		PSOCal();

	}




//DE Code
		for(i=0;i<D;i++)
		{
		printf("Dim%d\t\t",(i+1));
		}
		printf("\n");
	   	for (k=0; k < NO; k++)
		{
			printf("\nObjective Function NO::%d\n",(k+1));
			popul[k] = (int *)malloc(NP*(D+1)*sizeof(int));
			    //  if (popul[k][i] == NULL) perror("malloc");
			next[k] = (int *)malloc(NP*(D+1)*sizeof(int));
			      if (next[k] == NULL) perror("malloc");
			for (i=0; i < NP; i++)
			{
			      popul[k][i] = (int *)malloc((D+1)*sizeof(int));
			      if (popul[k][i] == NULL) perror("malloc");
				randm(arr,D);
			      for (j=0; j < D; j++)
				{
					popul[k][i][j] = arr[j];
					printf("%d\t",popul[k][i][j]);
				}
				if(k==0)
				{
			     	popul[k][i][D] = Distance(popul[k][i],D);
				printf("%d\t",popul[k][i][D]);
				printf("\n");
				}
				if(k==1)
				{
			     	popul[k][i][D] = NumberOfV(popul[k][i],D);
				//popul[k][i][D] = func(popul[k][i]);
				printf("%d\t",popul[k][i][D]);
				printf("\n");
				}
				if(k==2)
				{
			     	popul[k][i][D] = Makespan(popul[k][i],D);
				//popul[k][i][D] = func(popul[k][i]);
				printf("%d\t",popul[k][i][D]);
				printf("\n");
				}
				if(k==3)
				{
			     	popul[k][i][D] = WaitingTime(popul[k][i],D);
				//popul[k][i][D] = func(popul[k][i]);
				printf("%d\t",popul[k][i][D]);
				printf("\n");
				}
				if(k==4)
				{
			     	popul[k][i][D] =TotalDelay(popul[k][i],D);
				//popul[k][i][D] = func(popul[k][i]);
				printf("%d\t",popul[k][i][D]);
				printf("\n");
				}
			//printf("\n");
			     	numofFE++;

			     next[k][i] = (double *)malloc((D+1)*sizeof(double));
			      if (next[k][i] == NULL) perror("malloc");
			}//i
		}//k
//step 2,3


	copyLastPopPSO();
	for (j=0;j<Gmax;j++)
	{
		vectorCal();
	}


	endtime = clock();
   	totaltime = (double)(endtime - starttime);
//find the best solution
	findBestSol();
//Output
	display();
	freeMem();


}//end main

int copyLastPopPSO()
	{
	int k,i,j;
	for(k=0;k<NO;k++)
		{
		for(i=0;i<NP;i++)
			{
			for(j=0;j<D;j++)
				{
				popul[k][i][j]=populPSO[k][i][j];
				}
			}
		}
	}
int Distance(int *X,int M)
	{
	double xx=0.0;
	double dist1;
	int y11;
		int i;
		totalcost=0.0;
		for (i=1; i<M; i++)
	  	 {
			dist1=0;
			city1=X[i-1];
			city2=X[i];
			x1=mat[city1][1];
 			y11=mat[city1][2];
  			x2=mat[city2][1];
			y2=mat[city2][2];

			x=x1-x2;
			y=y11-y2;
			xx=x*x;
			yy=y*y;
			xy=xx+yy;
			dist1=sqrt(xy);
			totalcost=totalcost+dist1;
		 }

		return totalcost;
	}

//NumberOfV
int NumberOfV(int *X,int D)
	{
		double xx=0.0;
		int i,j;
		double Vcapacity=200;
		double vcapcitypick=0;
		 NumberOfVechiles=0;
		double pick;
		
		int cn=0;
		for (i=0; i<D; i++)
	  	 {

			city1=X[i];
			x1=mat[city1][3];
			pick=mat[city1][4];
			
			if(x1>Vcapacity)
			{
			citysToRoute[NumberOfVechiles]=cn;
			cn=0;
			Vcapacity=200;
			vcapcitypick=0;
			NumberOfVechiles++;
			Vcapacity=Vcapacity-x1;
			citysToRoute[NumberOfVechiles]=cn;
			}
			else
			{
 			Vcapacity=(Vcapacity-x1);
			vcapcitypick+=pick;
			if((Vcapacity+vcapcitypick)>200)
			{
			
			NumberOfVechiles++;
			cn=0;
			 
			vcapcitypick=0;
			Vcapacity=Vcapacity-x1;
			}
			routes[NumberOfVechiles][cn]=city1;
			cn++;
			}
			
		 }
		printf("==========%d",NumberOfVechiles);
		 for (i=0; i<NumberOfVechiles; i++)
       			 {
				VrootCapacity[i]=Distance(routes[i],citysToRoute[i]);
			    printf("root time=%d",VrootCapacity[i]);
			  
			}
		for (i=0; i<NumberOfVechiles; i++)
	  	 {
			printf("\nVechile no::%d",i);
			for (j=0; j<=citysToRoute[i]; j++)
	  		 {
				printf("\t %d",routes[i][j]);
                         }
		}
		printf("\n");

		return NumberOfVechiles;
	}
int Makespan(int *X,int D)
    {
        int i,j;
        int maxtime=0;
        printf("\n at makespan");
       for (i=0; i<NumberOfVechiles; i++)
        {

            printf("\n root time=%d",VrootCapacity[i]);
            if(VrootCapacity[i]>maxtime)
            {
                maxtime=VrootCapacity[i];
            }
	}
        printf("\n makespan= %d",maxtime);
       return maxtime;
    }
int WaitingTime(int *X,int D)
    {
         int i,XX=0;
	AA=0;
	double nn;
	double lc[1000],ac[1000],wc[1000],sc=90,dt[1000];
     	for(i=1;i<D;i++)
        {
	
		double xx=0.0;
		double dist1;
		int y11;
			dist1=0;
			city1=X[i-1];
			city2=X[i];
			x1=mat[city1][1];
 			y11=mat[city1][2];
  			x2=mat[city2][1];
			y2=mat[city2][2];

			x=x1-x2;
			y=y11-y2;
			xx=x*x;
			yy=y*y;
			xy=xx+yy;
			dist1=sqrt(xy);
	ac[X[i]]=lc[X[i-1]]+dist1;//2
	if(ac[X[i]]>=bc[X[i]]) //3
	{
	wc[X[i]]=0;
	}
	else
	{
	wc[X[i]]=bc[X[i]]-ac[X[i]];
	
	}
	
	

		
		
	lc[X[i]]=ac[X[i]]+wc[X[i]]+sc;	
	
	if(ac[X[i]]<=ec[X[i]]) //3
	{
	dt[X[i]]=0;
	}
	else
	{
	AA++;
	dt[X[i]]=ac[X[i]]-ec[X[i]];
	}
		
	
	
	
        XX=XX+wc[X[i]];
        }
	//printf("NN=%d",nn);
        return fabsf(XX);
    }

int TotalDelay(int *X,int D)
    {
	return fabsf(AA*10);

      }



int findPbest(int j)
	{
		int i,k;
		printf("Pbest is Called");
		if(j==0)
			{
			for (k=0; k < NOPSO; k++)	
	     		 {
				for (i=0;i<NPPSO;i++)
					{
					pbest[k][i]=populPSO[k][i][DPSO];
					}
			}
			}
		else
			{
			for(k=0;k<NOPSO;k++)
			{
			for (i=0;i<NPPSO;i++)
				{
				if(populPSO[k][i][DPSO]<pbest[k][i])
					{
					pbest[k][i]=populPSO[k][i][DPSO];
					}
				}
			}
			}
	}
int findGbest()
	{
		int i,k;
		printf("gbest is Called");
		for(k=0;k<NOPSO;k++)
		{
		gbest[k]=5000000000;
		for (i=0;i<NPPSO;i++)
			{
			if(pbest[k][i] < gbest[k])
				{
				gbest[k]=pbest[k][i];
				}
		
			}
		//printf("\ngbest is=%f",gbest);
		}
	}
int PSOCal()
	{
	int i,j,k,l,a,p,loc;
		
		int Uvalue[10];
		countPSO++;
		void swapUI(int *a,int *b)
		{
		int temp=*a;
		*a=*b;
		*b=temp;
		}
		
		printf("\n==================================================================");
		printf("\n\nTrial vectors U for iteration %d\n ",countPSO);
		for(k=0;k<DPSO;k++)
		{
			printf("U%d\t",(k+1));
			if(k==(DPSO-1))
			{
			printf("Fitness Value");	
			}
		}
		
	
		for (k=0; k < NOPSO; k++)	/* Going through whole population	*/
	     	 {
		printf("\n NOPSObjective Function NOPSO::%d\n",(k+1));
		for (i=0; i < NPPSO; i++)	/* Going through whole population	*/
		     	{
			r1 = URAND;
			r2 = URAND;
			r3 = URAND;
			
			printf("\n");
		for(a=0;a<DPSO;a++)
		{
		arrUI[a]=populPSO[k][i][a];
	//printf("%d\t",arrUI[a]);
		}
		printf("\n");
		
		sum=0;
		for (j=0; j < DPSO; j++)
		 {
			velocity[k][i][j] = W*velocity[k][i][j]+c1*r1*(pbest[k][i]-populPSO[k][i][j])+c2*r2*(gbest[k]-populPSO[k][i][j]);
			//printf("\t\tVelocity=%f\t",velocity[i][j]);
		
			UPSO[j]=populPSO[k][i][j]+velocity[k][i][j];
						if(UPSO[j]>Xl &&UPSO[j]<Xu)
						{
						UPSO[j]=UPSO[j];
						}
						else
						{
						UPSO[j] = populPSO[k][i][j];			
						}

			/*if(UPSO[j]>0 && UPSO[j]<=DPSO)
			{


				for(a=0;a<DPSO;a++)
				{

					if(arrUI[a]==UPSO[j])
					{
					loc=UPSO[j];
					swapUI(&arrUI[j],&arrUI[a]);

					}


				}


				//at last replace the U[j] by arrUI

			}
			else
			{
			j=j-1;
			r1 = (int)(NPPSO*URAND);
			 r2 = (int)(NPPSO*URAND);
			r3 = (int)(NPPSO*URAND);

			}
			if(j==(DPSO-1))
				{
					for(p=0;p<DPSO;p++)
					{
					UPSO[p]=arrUI[p];
					}
				}*/


			printf("%d\t",UPSO[j]);
		 }
				if(k==0)
				{
				UPSO[DPSO]=Distance(populPSO[k][i],DPSO);
				printf("%d",UPSO[DPSO]);
				Uvalue[k]=UPSO[DPSO];
				}
				if(k==1)
				{
				//U[D]=func(popul[k][i]);
				UPSO[DPSO]=NumberOfV(populPSO[k][i],DPSO);
				printf("%d",UPSO[DPSO]);
				Uvalue[k]=UPSO[DPSO];
				}
				if(k==2)
				{
				//U[D]=func(popul[k][i]);
				UPSO[DPSO]=Makespan(populPSO[k][i],DPSO);
				printf("%d",UPSO[DPSO]);
				Uvalue[k]=UPSO[DPSO];
				}
				if(k==3)
				{
				//U[D]=func(popul[k][i]);
				UPSO[DPSO]=WaitingTime(populPSO[k][i],DPSO);
				printf("%d",UPSO[DPSO]);
				Uvalue[k]=UPSO[DPSO];
				}
				if(k==4)
				{
				//U[D]=func(popul[k][i]);
				UPSO[DPSO]=TotalDelay(populPSO[k][i],DPSO);
				printf("%d",UPSO[DPSO]);
				Uvalue[k]=950;
				}

			sum=0;
				for(l=0;l<NOPSO;l++)
						{
						sum=sum+(weightCector[l]*Uvalue[l]);
						}
				printf("\nWeighted sum is=%f",sum);
				

		 numofFE++;
	 /* Comparing the trial vector 'U' and the old individual
		'next[i]' and selecting better one to continue in the
		   next population.	*/
		 if (UPSO[DPSO] <= populPSO[k][i][DPSO])
					 {
					    iptrPSO = UPSO;
					    UPSO = nextPSO[k][i];
					    nextPSO[k][i] = iptrPSO;
					}
					 else
					 {
						   for (j=0; j <= DPSO; j++)
							  nextPSO[k][i][j] = populPSO[k][i][j];
					}
				}//I
			
		       /* Pointers of old and new populations are swapped	*/
			}//K
			      ptrPSO = populPSO;
			      populPSO = nextPSO;
			      nextPSO = ptrPSO;
			
	return 0;

	}//end of PSOCall

int vectorCal()
	{	
		int i,j,k,l,a,p,loc;
		count++;
		int Uvalue[10];
		//swap U[i]
		void swapUI(int *a,int *b)
		{
		int temp=*a;
		*a=*b;
		*b=temp;
		}
		printf("\nTrial vectors U for iteration %d\n ",count);
		for(l=0;l<D;l++)
		{
			printf("U%d\t\t",(l+1));
		}
		printf("\n");
		for (k=0; k < NO; k++)	/* Going through whole population	*/
	     	 {
		printf("\nObjective Function NO::%d\n",(k+1));
			for (i=0; i < NP; i++)	/* Going through whole population	*/
		     	 {
				    r1 = (int)(NP*URAND);
			 r2 = (int)(NP*URAND);
			r3 = (int)(NP*URAND);
		
		 jrand = (int)(D*URAND);
	
		
		for(a=0;a<D;a++)
		{
		arrUI[a]=popul[k][i][a];
	//printf("%d\t",arrUI[a]);
		}
		printf("\n");

		sum=0;
		//mutation
		 for (j=0; j < D; j++)
		 {
			U[j]=0;
		   if (URAND <=CR || j == jrand)
		   {
		       U[j] = popul[k][r3][j] + F*(popul[k][r1][j] - popul[k][r2][j]);
			//printf("(%d)\t",U[j]);
			if(U[j]>0 && U[j]<=D)
			{


				for(a=0;a<D;a++)
				{

					if(arrUI[a]==U[j])
					{
					loc=U[j];
					swapUI(&arrUI[j],&arrUI[a]);

					}


				}


				//at last replace the U[j] by arrUI

			}
			else
			{
			j=j-1;
			r1 = (int)(NP*URAND);
			 r2 = (int)(NP*URAND);
			r3 = (int)(NP*URAND);

			}


		 }else
			{
			j=j-1;
			r1 = (int)(NP*URAND);
			 r2 = (int)(NP*URAND);
			r3 = (int)(NP*URAND);

			}
		   /* else
			{
	      		  U[j] = popul[i][j];

				for(a=0;a<D;a++)
				{
					if(arrUI[a]==U[j])
					{
					//swap
					swapUI(&arrUI[a],&arrUI[j]);
					}
				}

			}
			*/
			if(j==(D-1))
				{
					for(p=0;p<D;p++)
					{
					U[p]=arrUI[p];
					}
				}
			//printf("\n\n\n");

			//printf("%d\t",U[j]);
		}
		//U[D]=0;
		printf("\n");
		for(a=0;a<D;a++)
		{
		//arrUI[a]=popul[i][a];
		printf("%d\t",U[a]);
		}
				if(k==0)
				{
				U[D]=Distance(U,D);
				printf("%d",U[D]);
				Uvalue[k]=U[D];
				}
				if(k==1)
				{
				//U[D]=func(popul[k][i]);
				U[D]=NumberOfV(U,D);
				printf("%d",U[D]);
				Uvalue[k]=U[D];
				}
				if(k==2)
				{
				//U[D]=func(popul[k][i]);
				U[D]=Makespan(U,D);
				printf("%d",U[D]);
				Uvalue[k]=U[D];
				}
				if(k==3)
				{
				//U[D]=func(popul[k][i]);
				U[D]=WaitingTime(U,D);
				printf("%d",U[D]);
				Uvalue[k]=U[D];
				}
				if(k==4)
				{
				//U[D]=func(popul[k][i]);
				U[D]=TotalDelay(U,D);
				printf("%d",U[D]);
				Uvalue[k]=U[D];
				}
				sum=0;
				for(l=0;l<NO;l++)
						{
						sum=sum+(weightCector[l]*Uvalue[l]);
						}
				printf("\nWeighted sum is=%f",sum);
				
				
				
				 numofFE++;
					 /* Comparing the trial vector 'U' and the old individual
						'next[i]' and selecting better one to continue in the
							   next population.	*/
					 if (U[D] <= popul[k][i][D])
					 {
					    iptr = U;
					    U = next[k][i];
					    next[k][i] = iptr;
					}
					 else
					 {
						   for (j=0; j <= D; j++)
							  next[k][i][j] = popul[k][i][j];
					}
				}//I
			
		       /* Pointers of old and new populations are swapped	*/
			}//K
			      ptr = popul;
			      popul = next;
			      next = ptr;
			
	return 0;

   	}

int findBestSol()
	{
	int i,k;
		for (k=0; k < NO; k++)
	   	 {
	     	
		 for (i=0; i < NP; i++)
	   	 {
	     		 if (popul[k][i][D] < min_value)
	      		{
			 min_value = popul[k][i][D];
		 	index[k] = i;
	    		}
  		 }
		}
	}
	int display()
	{
	/* Printing out information about optimization process for the user	*/
		int i,k,m;
		printf("\n\nExecution time:::: %.3f s\n", (totaltime / (double)CLOCKS_PER_SEC));
		//printf("Number of Objective Function Evaluations: %d\n", numofFE);
		//printf("Solution:\nValues of variables:\n ");
		 for (k=0; k < NO; k++)
		{
		m=index[k];
		printf("\n---------------------------Objective Function Number= %d ----------------------------------------------\n",k);
		 for (i=0; i < D; i++)
			{
		      printf("%d \t ", popul[k][m][i]);
			}
			printf("\nObjective function value %d=",k);



		  	 printf("%d\n", popul[k][m][D]);
		}
		m=index[0];
		fprintf(fod,"\nDistance =%d\n",popul[0][m][D]);
		m=index[1];
		fprintf(fod,"\nNumber of Vehicles =%d\n",popul[1][m][D]);
		m=index[2];
		fprintf(fod,"\nMakespan =%d\n",popul[2][m][D]);
		m=index[3];
		fprintf(fod,"\nwaiting Time =%d\n",popul[3][m][D]);
		m=index[4];
		fprintf(fod,"\nDelay Time =%d\n",popul[4][m][D]);

		fprintf(fod,"\nWeighted Sum Approach=%f\n",sum);
		fprintf(fod,"\nTotal Time=%f\n",(totaltime / (double)CLOCKS_PER_SEC));
		for (i=0; i<NumberOfVechiles; i++)
	  		 {
			fprintf(fod,"\nVechile no::%d\t=",i);
			fprintf(fod,"0");
				for (j=0; j<=citysToRoute[i]; j++)
		  		 {
					fprintf(fod,"\t %d",routes[i][j]);
		                 }
			fprintf(fod,"\t0");
			}

	}
	int freeMem()
	{
		int i,k;
		 for (k=0; k < NOPSO; k++)
		 {
		 for (i=0; i < NPPSO; i++)
		 {
		 	free(populPSO[k][i]);
		      	free(nextPSO[k][i]);
			free(popul[k][i]);
		      	free(next[k][i]);
		  }
		}
		free(populPSO);
		free(nextPSO);
		free(UPSO);
		free(popul);
		free(next);
		free(U);
	}
