#include<stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

int main(void)
{
    MPI_Init(NULL, NULL);

    int i,j,k;
    const int n=8192;
    int map[n];

    const int total_num_tests=5;

    double sum=0.0;
    double* A;

    double b[n],c[n], x[n];
    A= malloc(n*n*sizeof(double));


    double* A_copia;
    double* b_copia;


    double result_tests[total_num_tests];
    double temp;

    int rank;
    int comm_sz;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    double start, finish, loc_elapsed, elapsed;
    MPI_Comm comm;

    comm=MPI_COMM_WORLD;



   for(int num_test=0; num_test <total_num_tests; num_test++) {
    if (rank==0)
    {

      A_copia=  malloc(n*n*sizeof(double));
      b_copia= malloc(n*sizeof(double));

              for (i = 0; i < n; i++){
                 sum=0.0;
                 b[i] = (10.0 * (double) rand())/((double) RAND_MAX);
                 b_copia[i]= b[i];
                 for (j = 0; j < n; j++){

                   A[i*n + j] = (10 * ((double) rand()))/((double) RAND_MAX);
                   A_copia[i*n + j]= A[i*n + j];
                   sum =sum +A[i*n+j];
                   if (j == n-1) {
                      A[i*n+i]=sum ;
                      A_copia[i*n+i]=A[i*n+i];
                   }

                 }
              }

    }




    MPI_Bcast (A,n*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast (b,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
   start = MPI_Wtime();
    for(i=0; i<n; i++)
    {
        map[i]= i % comm_sz;
    }

    for(k=0;k<n;k++)
    {
        MPI_Bcast (&A[k*n + k],n-k,MPI_DOUBLE,map[k],MPI_COMM_WORLD);
        MPI_Bcast (&b[k],1,MPI_DOUBLE,map[k],MPI_COMM_WORLD);
        for(i= k+1; i<n; i++)
        {
            if(map[i] == rank)
            {
                c[i]=A[i*n+k]/A[k*n+k];
                for(j=k;j<n;j++)
                {
                    A[i*n+j]=A[i*n+j]-( c[i]*A[k*n+j] );
                }
                b[i]=b[i]-( c[i]*b[k] );
            }
        }

    }






    if (rank==0)
    {
        //Si se puede incluir dentro del lazo pero eso solo una operacion
        // y queda más claro dejarlo fuera del lazo
        x[n-1]=b[n-1]/A[(n-1)*n + n-1];
        for(i=n-2;i>=0;i--)
        {
            sum=0;

            for(j=i+1;j<n;j++)
            {
                sum=sum+A[i*n+j]*x[j];
            }
            x[i]=(b[i]-sum)/A[i*n + i];
        }


    }

    finish = MPI_Wtime();
    loc_elapsed = finish-start;
    MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    result_tests[num_test]=elapsed;

    if (rank==0)
    {

        double res= 0.0;
        for(i=0;i<n;i++)
        {   temp=0.0;
            for(j=0;j<n;j++){
                temp=temp+A_copia[i*n + j]*x[j];

            }
            res+= fabs(b_copia[i]-temp);

        }
        printf("\nEl residuo  es: %.2lf\n", res);


        printf("tamano de la matriz = %d\n", n);
        printf("num procesos = %d\n", comm_sz);
        printf("Tiempo transcurrido experimento= %e(s)\n", elapsed);

    }


  }

  if(rank==0) {
  sum=0.0;
  for(i=0;i<total_num_tests;i++){
    sum=sum+result_tests[i];
  }

  printf("tamano de la matriz = %d\n", n);
  printf("num procesos = %d\n", comm_sz);
  printf("Tiempo transcurrido promedio = %e(s)\n", sum / total_num_tests);
  }

  MPI_Finalize();

  free(A);
  free(A_copia);
  free(b_copia);
  return 0;


}
