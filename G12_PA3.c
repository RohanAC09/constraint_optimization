#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

double Conjugate_direction(int p);
double objective_function(int p, int c, double X[c]);
double Bounding_Secant(int p, int c, int k, double X[c], double S[c][c], double R);
double objective(int p, int c, double y, int k, double X[c], double S[c][c], double R);
double df(int p, int c, double y, int k, double X[c], double S[c][c], double R);
double Gauss(int c, double A[c][c]);
double Constraints(int p, int c, double X[c]);
double Penalty(int p, int c, double X[c], double R);

double objective_function(int p, int c, double X[c])		//Objective Function
{
double f=0;
switch(p)
{	case 1:														//Problem 1
		f  = pow( X[0] -10 ,3) + pow( X[1] -20 ,3);				// Main function
		return(f);		break;
	case 2: 																	//Problem 1
		f = pow(sin(2*M_PI*X[0]),3)*sin(2*M_PI*X[1])/(pow(X[0],3)*(X[0]+X[1])) ;	// Main function
		return(f);		break;
	case 3:												//Dixon-Price Function
		f = X[0] + X[1] + X[2] ;
		return(f);		break;
}}

double Constraints(int p, int c, double X[c])
{
double q[25];
int i;
switch(p)
{	case 1:														//Problem 1
		q[0] = (pow( X[0] -5 ,2) + pow( X[1] -5 ,2) -100)/150 ;		// Constraints
		q[1] = (82.81 - pow( X[0] -6 ,2) - pow( X[1] -5 ,2))/32.81 ;
		q[2] = (X[0] -13 )/7;
		q[3] = (20 - X[0])/7 ;
		q[4] = (4 - X[1])/4 ;
		q[5] = (X[1])/4;
		q[6] = 0;
		for(i=0;i<6;i++)
			if(q[i]>0)
				q[i]=0;
		for(i=0;i<6;i++)
			q[6] = q[6] + pow(q[i],2) ;
		return(q[6]);		break;
	case 2: 																	//Problem 1
		q[0] = (X[1] - pow( X[0],2) -1)/9 ;						// Constraints
		q[1] = (X[0] - pow( X[1] -4 ,2) - 1)/9 ;
		q[2] = (10 - X[0])/10 ;
		q[3] = (10 - X[1])/10 ;
		q[4] = (X[0])/10;
		q[5] = (X[1])/10;
		q[6] = 0;
		for(i=0;i<6;i++)
			if(q[i]>0)
				q[i]=0;
		for(i=0;i<6;i++)
			q[6] = q[6] + pow(q[i],2) ;
		return(q[6]);		break;
	case 3:	
		q[0] = (1 - 0.0025*(X[3]+X[5]))*1000 ;						// Constraints
		q[1] = (1 - 0.0025*(X[6]+X[4]-X[3]))*300 ;
		q[2] = (1 - 0.01*(X[7]-X[5]))*100 ;
		q[3] = (X[0]*X[5] - 100*X[0] - 833.33252*X[3] + 83333.333)/10000 ;
		q[4] = (X[1]*X[6]-X[1]*X[3]+1250*X[3]-1250*X[4])/11000 ;
		q[5] = (X[2]*X[7]-X[2]*X[4]+2500*X[4]-1250000)/11000 ;
		q[6] = (10000 - X[0])/9.9 ;
		q[7] = (10000 - X[1])/9 ;
		q[8] = (10000 - X[2])/9 ;
		q[9] = (X[0] -100)/9.9 ;
		q[10] = (X[1] -1000)/9 ;
		q[11] = (X[2] -1000)/9 ;
		q[22]=0;
		for(i=3;i<8;i++)
			q[i+9] = (1000 - X[i]) ;
		for(i=3;i<8;i++)
			q[i+14] = (X[i] -10) ;
		for(i=0;i<22;i++)
			if(q[i]>0)
				q[i]=0;
		for(i=0;i<22;i++)
			q[22] = q[22] + pow(q[i],2) ;
		return(q[22]);		break;
}
}

double Penalty(int p, int c, double X[c], double R)
{
double f;
switch(p)
{	case 1:														//Problem 1
		f = objective_function(p,c,X) + R* Constraints(p,c,X);
		return(f);		break;
	case 2: 																	//Problem 1
		f = -objective_function(p,c,X) + R* Constraints(p,c,X);
		return(f);		break;
	case 3:	
		f = objective_function(p,c,X) + R* Constraints(p,c,X);
		return(f);		break;
}
}

main()													// Main Function body
{
	int p;												// p = problem number from the given set
	printf("\n**********************************\nPenalty Function Method\n");
	printf("\n**********************************\n");
	printf("\nPlease enter problem No. = ");
	scanf("%d",&p);
	
	Conjugate_direction(p);
	return 0;						
}

double Conjugate_direction(int p)
{
	FILE *f;
	switch(p)
	{
	case 1:	f=fopen("In_P1.txt","r"); break;
	case 2:	f=fopen("In_P2.txt","r"); break;
	case 3:	f=fopen("In_P3.txt","r"); break;
	}
	double a,b,E1,E2,R;										//a= UL, b= LL, E1 = Accuracy
	int c,i,j,M,k=1,N,Ns=0,l,u,q,r=10;						//M = Max Iterations & c = Number of variables (Input)
	fscanf(f,"No of var=%d\tLL=%lf\tUL=%lf\nAccuracy1=%lf\tAccuracy2=%lf\tMax_Itr=%d\n\n",&c,&a,&b,&E1,&E2,&M);
	double S[c][c],X[c],alpha,normDir;					//S = directions & X = coordinates
	double X1[c],S1[c],S2[c][c],z,G,fun,f1,f2;				//(in current itr) X1 = 1st var set, S1 = 1st dir set
	l=a,u=b;
	fclose(f);
	srand(time(0));
	switch(p)
	{
		case 1: f=fopen("Out_P1.txt","w"); break;
		case 2: f=fopen("Out_P2.txt","w"); break;
		case 3: f=fopen("Out_P3.txt","w"); break;
	}
/**/	for(q=0;q<10;q++)
	{
	for(i=0;i<c;i++)
	{
		for(j=0;j<c;j++)								//Initializing unit directions
		{
			if(i==j)
			S[i][j]=1;
			else
			S[i][j]=0;
		}
		X[i]=(rand()%( u - l +1)) + l;					//Storing Initial guess in 1st col of X
/**/	}k=1,Ns=0,R=0.1;
	printf("\nInitial Guess :\n\t");					//printing initial guess
	for(i=0;i<c;i++)
	printf("X_%d = %lf\t\t",i+1,X[i]);
	fprintf(f,"Itr\t");									// output file data
	for(i=0;i<c;i++)
	fprintf(f,"X_%d\t\t",i+1);
	fprintf(f,"Fun Value\tSearches\n\n 0\t");
	for(i=0;i<c;i++)
	fprintf(f,"%lf\t",X[i]);
	fun=objective_function(p,c,X);						//function evaluation for initial guess
/**/	f1=fun;
	fprintf(f,"%lf\t%d\t",fun,Ns);
	normDir=100; N=c*c ;								// norm of vector d
	while(k<(M+1))											//Starting iterative loop
	{
		fprintf(f,"\n %d\t",k);
		for(i=0;i<c;i++)
		{
			alpha = Bounding_Secant(p,c,i,X,S,R);			// N unidirectional search
			for(j=0;j<c;j++)
				X[j] = X[j] + alpha*S[j][i] ;
			if(i==0)									// Storing X1 for 'd' vector (d=Xn-X1)
			for(j=0;j<c;j++)
				X1[j] = X[j];
			Ns++;										// Count of unidirectional searches
		}
		alpha = Bounding_Secant(p,c,0,X,S,R);				// Last S1 search
		for(j=0;j<c;j++)
			X[j] = X[j] + alpha*S[j][0] ;
		for(i=0;i<c;i++)
		fprintf(f,"%lf\t",X[i]);
		fun=objective_function(p,c,X); Ns++;
		fprintf(f,"%.8lf\t%d\t",fun,Ns);
		f2=f1-fun;										// Function Accuracy condition
/**/	if(fabs(f2)<(E2/10000))
		break;
		f1=fun;
		for(j=0;j<c;j++)								// 'd' vector
			S1[j] = X[j] - X1[j] ;
		z=0;
		for(i=0;i<c;i++)
			z = z + pow(S1[i],2) ;
		normDir = sqrt(z);
		for(j=c-1;j>=0;j--)
			for(i=0;i<c;i++)
				S[i][j] = S[i][j-1] ;					// Changing directions (Shifting S3=S2, S2=S1 ...)
		for(i=0;i<c;i++)
			S[i][0] = (S1[i])/(normDir) ;				// updated 1st search direction
		k++;
		for(j=0;j<c;j++)
			for(i=0;i<c;i++)
				S2[i][j] = S[i][j] ;
		G=Gauss(c,S2);									// Checking 'Linear Independency'
/**/	if(normDir<E1 || fabs(f2)<E1 || G==0)					// Checking norm condition
		for(i=0;i<c;i++) 
		{
			for(j=0;j<c;j++)							//Re-Initializing unit directions
			{
				if(i==j)
				S[i][j]=1;
				else
				S[i][j]=0;
			}}	R=R*r;
	}
	printf("\n\nFinal Solution is \n\n\t");
	for(i=0;i<c;i++)
		printf("X_%d = %lf\t\t",i+1,X[i]);
	printf("\n\n");
	fprintf(f,"\n\n");
	}
	fclose(f);
}

double Bounding_Secant(int p, int c, int k, double X[c], double S[c][c], double R)		/* Bounding Phase */
{
	double delta=0.0000001,alpha;						// Icrement value & limits
	double x0=0, x1, x2, f0, f1, f2,te;
	int i=2;
	x1 = (x0 - delta) ; x2 = (x0 + delta) ;			// New points
	f0 = objective(p,c,x0,k,X,S,R);					// Calculate objective_function
	f1 = objective(p,c,x1,k,X,S,R);
	f2 = objective(p,c,x2,k,X,S,R);
	if(f1 <= f0)									// Deciding the direction and sign of Delta
		if(f0 <= f2)
			delta= -delta;
	else
		x0 = x1, x1 = x2, f1 = f2;					// Declaring first value	
	
	x2 = x1 + 2*delta;								/* Step 2 */
	f2 = objective(p,c,x2,k,X,S,R);
	do{
		x0 = x1;	x1 = x2;	x2 = x1 + pow(2,i)*delta;		/*If not terminated,then update values*/
		f1 = f2;	f2 = objective(p,c,x2,k,X,S,R);
		i++;
	}while(f1 > f2);								//while number lies within limits
	if(x2<x0){
		te=x2,	x2=x0,	x0=te;							// Interchanging values
	}
	double z,f3,f4=0,f5,e=0.001;
    f1=df(p,c,x0,k,X,S,R),	f2=df(p,c,x2,k,X,S,R);		//Gradient
	do{
        if(fabs(f1-f2)<e)										//Checking gradients
    		{z=(x0+x2)/2;break;}
        if(f1>0 || f2<0)								//condition of unapplicable method
        	{z=(x0+x2)/2;break;}
		z = x2 - ((f2 * (x2-x0)) / (f2-f1));			//Calculating next term
        f3=df(p,c,z,k,X,S,R);
		if(f3<0)
		x0=z,	f1=f3;
		else
		x2=z,	f2=f3;
		f5=f3-f4;
		if(fabs(f5)<e)
			break;
		f4=f3;
    } while(fabs(f3)>e);								//End of iteration's loop
	return(z);
}

double df(int p, int c, double y, int k, double X[c], double S[c][c], double R)	//Derivative using central difference method
{
	double m=0.000001;
	double dy = (objective(p,c,(y+m),k,X,S,R)-objective(p,c,(y-m),k,X,S,R))/(2*m);
    return(dy);	
}

double objective(int p, int c, double y, int k, double X[c], double S[c][c], double R)	// X + (Alpha)*S
{
	double f,A[c];
	int i;
	for(i=0;i<c;i++)
		A[i] = X[i] + y*S[i][k];
	f = Penalty(p,c,A,R);
	return(f);
}

double Gauss(int c, double A[c][c])						// Linear independency by Gauss matrix method
{
	int i,j,k,l,z;
	double temp,a;
	for(j=0;j<c;j++)						//Forward Elimination
	{
		temp= fabs(A[j][j]);				//Pivot value
		z=j;
		for(l=j+1;l<c;l++)
			if(fabs(A[l][j])>temp)			//larger (non-zero) number
			{
				z=l;						//Row number
				temp=A[l][j];
			}
		if(z!=j)
			for(k=0;k<c;k++)				//interchanging rows
			{
				temp=A[j][k];
				A[j][k]=A[z][k];
				A[z][k]=temp;
			}
		if(A[j][j]==0)						// Zero column condition
		{
			return 0; break;
		}
		for(i=0;i<c;i++)
		{
			if(i>j)
			a=A[i][j]/A[j][j];
			for(k=0;k<c;k++)
				A[i][k]=A[i][k]-a*A[j][k];	//Row transformation
		}
	}
	if(fabs(A[c-1][c-1])<0.0001)
	return(0);								// Dependent
	else
	return(1);								// Independent
}
