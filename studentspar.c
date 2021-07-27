/*
Grupo - 01
Alunos:
João Pedro Fidelis Belluzzo - 10716661
Leonardo Cerce Guioto - 10716640
Luis Fernando Costa de Oliveira - 10716532
Rodrigo Augusto Valeretto - 10684792
Thiago Daniel Cagnoni de Pauli - 10716629
*/

/*
Compilação e execução: make -f Makefile.par < input.txt
Compilação: mpicc studentspar.c -o studentspar -fopenmp -lm
Execução: mpirun -np 6 --hostfile hostfile.txt --map-by node ./studentspar < input.txt
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <mpi.h>
#include <limits.h>
#include <locale.h>

#define T 8

typedef struct{
    int * min;
    int * max;
    double * mean;
    double * sd;
    double * median;
}analises;

int find_max(int * vec, int tam){
    int maximo = INT_MIN;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(max : maximo)
    for(int i=0; i<tam; i++){
        if(vec[i] > maximo){
            maximo = vec[i];
        }
    }

    return maximo;
}

int find_max_index(double * vec, int tam){
    double maximo = LONG_MIN;
    int max_index;

    #pragma omp parallel for num_threads(T) shared(vec, max_index)
    for(int i=0; i<tam; i++){
        #pragma omp critical
        {
            if(vec[i] > maximo){
                maximo = vec[i];
                max_index = i;
            }
        }
            
        
    }

    return max_index;
}

int find_min(int * vec, int tam){
    int minimo = INT_MAX;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(min : minimo)
    for(int i=0; i<tam; i++){
        minimo = (minimo > vec[i]) ? vec[i] : minimo;
    }

    return minimo;
}

double find_city_mean(int * vec, int tam){
    double mean = 0;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(+ : mean)
    for(int i=0; i<tam; i++){
        mean += vec[i];
    }

    mean /= tam;

    return mean;
}

double find_mean(double * vec, int tam){
    double mean = 0;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(+ : mean)
    for(int i=0; i<tam; i++){
        mean += vec[i];
    }

    mean /= tam;

    return mean;
}

double find_sd(int * vec, int tam, double mean){
    double sd = 0;

    #pragma omp parallel for num_threads(T/2) shared(vec) reduction(+ : sd)
    for(int i=0; i<tam; i++){
        sd += pow(((double) vec[i] - mean),2);
    }

    sd /= tam;
    sd = sqrt(sd);

    return sd;
}

//Função de comparação utilizada no quicksort
int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int main(int argc, char * argv[]){

    int provided;
    setlocale(LC_ALL, "Portuguese");

    // Inicia o MPI indicando a utilização de threads com OMP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    // Se não estiver disponível a utilização de múltiplas threads, encerra o processo
    if(provided != MPI_THREAD_MULTIPLE){
        printf("Problema na criação das threads\n");
        return 1;
    }

    int seed; //Valores de entrada
    int best_mean_region_index, best_city_index_i, best_city_index_j;
    int *grades;
    double time;
    double best_mean_region, best_mean_city; 

    int root = 0, dst = 0, tag = 0, src;
    int my_rank, num_proc;
    int RCA[3];

    omp_set_nested(1);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    
    if(my_rank == 0){
        scanf("%d", &(RCA[0]));
        scanf("%d", &(RCA[1]));
        scanf("%d", &(RCA[2]));
        scanf("%d", &seed);
    }

    MPI_Bcast(RCA, 3, MPI_INT, root, MPI_COMM_WORLD);

    grades = (int*) malloc(RCA[0] * RCA[1] * RCA[2] * sizeof(int));

    if(my_rank == 0){
        //Inicia a geração das notas com base na seed da entrada
        srand(seed);

        for(int i=0; i < RCA[0] * RCA[1] * RCA[2]; i++){
            grades[i] = rand() % 101;
        }
    }

    MPI_Bcast(grades, RCA[0] * RCA[1] * RCA[2], MPI_INT, root, MPI_COMM_WORLD);

    time = omp_get_wtime();

    if(my_rank == 0){
        MPI_Request requests[15];
        MPI_Status status[15];
        int best_city, best_reg;

        analises analises_cidades;
        analises_cidades.min = (int *) malloc(RCA[0] * RCA[1] * sizeof(int));
        analises_cidades.max = (int *) malloc(RCA[0] * RCA[1] * sizeof(int));
        analises_cidades.mean = (double *) malloc(RCA[0] * RCA[1] * sizeof(double));
        analises_cidades.median = (double *) malloc(RCA[0] * RCA[1] * sizeof(double));
        analises_cidades.sd = (double *) malloc(RCA[0] * RCA[1] * sizeof(double));

        analises analises_regioes;
        analises_regioes.min = (int *) malloc(RCA[0] * sizeof(int));
        analises_regioes.max = (int *) malloc(RCA[0] * sizeof(int));
        analises_regioes.mean = (double *) malloc(RCA[0] * sizeof(double));
        analises_regioes.median = (double *) malloc(RCA[0] * sizeof(double));
        analises_regioes.sd = (double *) malloc(RCA[0] * sizeof(double));

        analises analise_brasil;
        analise_brasil.min = (int *) malloc(sizeof(int));
        analise_brasil.max = (int *) malloc(sizeof(int));
        analise_brasil.mean = (double *) malloc(sizeof(double));
        analise_brasil.median = (double *) malloc(sizeof(double));
        analise_brasil.sd = (double *) malloc(sizeof(double));

        src = 1; tag = 0;
        MPI_Irecv(analises_cidades.max, RCA[0] * RCA[1], MPI_INT, src, tag, MPI_COMM_WORLD, &(requests[0]));
        src = 2; tag = 0;
        MPI_Irecv(analises_cidades.min, RCA[0] * RCA[1], MPI_INT, src, tag, MPI_COMM_WORLD, &(requests[1]));
        src = 3; tag = 0;
        MPI_Irecv(analises_cidades.mean, RCA[0] * RCA[1], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &(requests[2]));
        src = 4; tag = 0;
        MPI_Irecv(analises_cidades.sd, RCA[0] * RCA[1], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &(requests[3]));
        src = 5; tag = 0;
        MPI_Irecv(analises_cidades.median, RCA[0] * RCA[1], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &(requests[4]));


        src = 1; tag = 1;
        MPI_Irecv(analises_regioes.max, RCA[0], MPI_INT, src, tag, MPI_COMM_WORLD, &(requests[5]));
        src = 2; tag = 1;
        MPI_Irecv(analises_regioes.min, RCA[0], MPI_INT, src, tag, MPI_COMM_WORLD, &(requests[6]));
        src = 3; tag = 1;
        MPI_Irecv(analises_regioes.mean, RCA[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &(requests[7]));
        src = 4; tag = 1;
        MPI_Irecv(analises_regioes.sd, RCA[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &(requests[8]));
        src = 5; tag = 1;
        MPI_Irecv(analises_regioes.median, RCA[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &(requests[9]));    


        src = 1; tag = 2;
        MPI_Irecv(analise_brasil.max, 1, MPI_INT, src, tag, MPI_COMM_WORLD, &(requests[10]));
        src = 2; tag = 2;
        MPI_Irecv(analise_brasil.min, 1, MPI_INT, src, tag, MPI_COMM_WORLD, &(requests[11]));
        src = 3; tag = 2;
        MPI_Irecv(analise_brasil.mean, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &(requests[12]));
        src = 4; tag = 2;
        MPI_Irecv(analise_brasil.sd, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &(requests[13]));
        src = 5; tag = 2;
        MPI_Irecv(analise_brasil.median, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &(requests[14]));       

        MPI_Wait(&(requests[0]), &(status[0])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[1]), &(status[1])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[2]), &(status[2])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[3]), &(status[3])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[4]), &(status[4])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[5]), &(status[5])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[6]), &(status[6])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[7]), &(status[7])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[8]), &(status[8])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[9]), &(status[9])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[10]), &(status[10])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[11]), &(status[11])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[12]), &(status[12])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[13]), &(status[13])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[14]), &(status[14])); // agora espera terminar a comm coletiva nao bloqueante

        best_city = find_max_index(analises_cidades.mean, RCA[0] * RCA[1]);
        best_reg = find_max_index(analises_regioes.mean, RCA[0]);

        time = omp_get_wtime() - time;

        for(int i=0; i<RCA[0]; i++){
            for(int j=0; j<RCA[1]; j++){
                printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2lf, média: %.2lf e DP: %.2lf\n", 
                i, j, analises_cidades.min[i*RCA[1] + j], analises_cidades.max[i*RCA[1] + j],
                analises_cidades.median[i*RCA[1] + j], analises_cidades.mean[i*RCA[1] + j],
                analises_cidades.sd[i*RCA[1] + j]);
            }
            printf("\n");
        }

        for(int i=0; i<RCA[0]; i++){
            printf("Reg %d: menor: %d, maior: %d, mediana: %.2lf, média: %.2lf e DP: %.2lf\n", 
                i, analises_regioes.min[i], analises_regioes.max[i],
                analises_regioes.median[i], analises_regioes.mean[i],
                analises_regioes.sd[i]);
        }

        printf("\n");

        printf("Brasil: menor: %d, maior: %d, mediana: %.2lf, média: %.2lf e DP: %.2lf\n", 
                *analise_brasil.min, *analise_brasil.max,
                *analise_brasil.median, *analise_brasil.mean,
                *analise_brasil.sd);

        printf("\n");

        printf("Melhor região: Região %d\n", best_reg);
        printf("Melhor cidade: Região %d, Cidade %d\n", best_city/RCA[1], best_city%RCA[1]);

    
        printf("\n");
        printf("Tempo de resposta sem considerar E/S, em segundos: %.8fs\n",time);
    }

    if(my_rank == 1){
        MPI_Request requests[3];
        MPI_Status status[3];
        dst = 0;

        int * max_cidades = (int*) malloc(RCA[0] * RCA[1] * sizeof(int));
        int * max_regioes = (int*) malloc(RCA[0] * sizeof(int));
        int max_brasil;

        #pragma omp parallel for num_threads(2) shared(max_cidades)
        for(int i=0; i<RCA[0] * RCA[1]; i++){
            max_cidades[i] = find_max(&(grades[i * RCA[2]]), RCA[2]);
            //printf("Cidade: %d | Max: %d\n",i,max_cidades[i]);
        }

        tag = 0;
        MPI_Isend(max_cidades, RCA[0] * RCA[1], MPI_INT, dst, tag, MPI_COMM_WORLD, &(requests[0]));

        #pragma omp parallel for num_threads(2) shared(max_regioes)
        for(int i=0; i<RCA[0]; i++){
            max_regioes[i] = find_max(&(max_cidades[i*RCA[1]]), RCA[1]);
            //printf("Região: %d | Max: %d\n",i,max_regioes[i]);
        }

        tag = 1;
        MPI_Isend(max_regioes, RCA[0], MPI_INT, dst, tag, MPI_COMM_WORLD, &(requests[1]));

        max_brasil = find_max(max_regioes, RCA[0]);
        //printf("Brasil | Max: %d\n",max_brasil);

        tag = 2;
        MPI_Isend(&max_brasil, 1, MPI_INT, dst, tag, MPI_COMM_WORLD, &(requests[2]));

        MPI_Wait(&(requests[0]), &(status[0])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[1]), &(status[1])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[2]), &(status[2])); // agora espera terminar a comm coletiva nao bloqueante

        free(max_cidades);
        free(max_regioes);
    }

    if(my_rank == 2){
        MPI_Request requests[3];
        MPI_Status status[3];

        dst = 0;
        int * min_cidades = (int*) malloc(RCA[0] * RCA[1] * sizeof(int));
        int * min_regioes = (int*) malloc(RCA[0] * sizeof(int));
        int min_brasil;

        #pragma omp parallel for num_threads(2) shared(min_cidades)
        for(int i=0; i<RCA[0] * RCA[1]; i++){
            min_cidades[i] = find_min(&(grades[i * RCA[2]]), RCA[2]);
            //printf("Cidade: %d | Min: %d\n",i,min_cidades[i]);
        }

        tag = 0;
        MPI_Isend(min_cidades, RCA[0] * RCA[1], MPI_INT, dst, tag, MPI_COMM_WORLD, &(requests[0]));

        #pragma omp parallel for num_threads(2) shared(min_regioes)
        for(int i=0; i<RCA[0]; i++){
            min_regioes[i] = find_min(&(min_cidades[i*RCA[1]]), RCA[1]);
            //printf("Região: %d | Min: %d\n",i,min_regioes[i]);
        }

        tag = 1;
        MPI_Isend(min_regioes, RCA[0], MPI_INT, dst, tag, MPI_COMM_WORLD, &(requests[1]));

        min_brasil = find_min(min_regioes, RCA[0]);
        //printf("Brasil | Min: %d\n",min_brasil);

        tag = 2;
        MPI_Isend(&min_brasil, 1, MPI_INT, dst, tag, MPI_COMM_WORLD, &(requests[2]));

        MPI_Wait(&(requests[0]), &(status[0])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[1]), &(status[1])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[2]), &(status[2])); // agora espera terminar a comm coletiva nao bloqueante

        free(min_cidades);
        free(min_regioes);

    }

    if(my_rank == 3){
        MPI_Request requests[5];
        MPI_Status status[5];

        double * mean_cidades = (double*) malloc(RCA[0] * RCA[1] * sizeof(double));
        double * mean_regioes = (double*) malloc(RCA[0] * sizeof(double));
        double mean_brasil;
        MPI_Request request;

        #pragma omp parallel for num_threads(2) shared(mean_cidades)
        for(int i=0; i<RCA[0] * RCA[1]; i++){
            mean_cidades[i] = find_city_mean(&(grades[i * RCA[2]]), RCA[2]);
            //printf("Cidade: %d | Mean: %.2lf\n",i,mean_cidades[i]);
        }

        dst = 4; tag = 1;
        MPI_Isend(mean_cidades, RCA[0] * RCA[1], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[0]));
        dst = 0; tag = 0;
        MPI_Isend(mean_cidades, RCA[0] * RCA[1], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[1]));


        #pragma omp parallel for num_threads(2) shared(mean_regioes)
        for(int i=0; i<RCA[0]; i++){
            mean_regioes[i] = find_mean(&(mean_cidades[i*RCA[1]]), RCA[1]);
            //printf("Região: %d | Mean: %.2lf\n",i,mean_regioes[i]);
        }

        dst = 4; tag = 2;
        MPI_Isend(mean_regioes, RCA[0], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[2]));
        dst = 0; tag = 1;
        MPI_Isend(mean_regioes, RCA[0], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[3]));

        mean_brasil = find_mean(mean_regioes, RCA[0]);
        //printf("Brasil | Mean: %.2lf\n",mean_brasil);

        dst = 4; tag = 3;
        MPI_Send(&mean_brasil, 1, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
        dst = 0; tag = 2;
        MPI_Isend(&mean_brasil, 1, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[4]));

        MPI_Wait(&(requests[0]), &(status[0])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[1]), &(status[1])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[2]), &(status[2])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[3]), &(status[3])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[4]), &(status[4])); // agora espera terminar a comm coletiva nao bloqueante

        free(mean_cidades);
        free(mean_regioes);

    }
    if(my_rank == 4){
        MPI_Request requests[3];
        MPI_Status statuss[3];

        double * mean_cidades = (double*) malloc(RCA[0] * RCA[1] * sizeof(double));
        double * sd_cidades = (double*) malloc(RCA[0] * RCA[1] * sizeof(double));
        double * mean_regioes = (double*) malloc(RCA[0] * sizeof(double));
        double * sd_regioes = (double*) malloc(RCA[0] * sizeof(double));
        double mean_brasil;
        double sd_brasil;
        
        MPI_Status status;

        src = 3; tag = 1;
        MPI_Recv(mean_cidades, RCA[0] * RCA[1], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

        #pragma omp parallel for num_threads(2) shared(sd_cidades)
        for(int i=0; i<RCA[0] * RCA[1]; i++){
            sd_cidades[i] = find_sd(&(grades[i * RCA[2]]), RCA[2], mean_cidades[i]);
            //printf("Cidade: %d | SD: %.2lf\n",i,sd_cidades[i]);
        }

        dst = 0; tag = 0;
        MPI_Isend(sd_cidades, RCA[0] * RCA[1], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[0]));

        src = 3; tag = 2;
        MPI_Recv(mean_regioes, RCA[0], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

        #pragma omp parallel for num_threads(2) shared(mean_regioes)
        for(int i=0; i<RCA[0]; i++){
            sd_regioes[i] = find_sd(&(grades[i*RCA[1]*RCA[2]]), RCA[1]*RCA[2], mean_regioes[i]);
            //printf("Região: %d | SD: %.2lf\n",i,sd_regioes[i]);
        }

        dst = 0; tag = 1;
        MPI_Isend(sd_regioes, RCA[0], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[1]));

        src = 3; tag = 3;
        MPI_Recv(&mean_brasil, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);

        sd_brasil = find_sd(grades, RCA[0] * RCA[1] * RCA[2], mean_brasil);
        //printf("Brasil | SD: %.2lf\n",sd_brasil);

        dst = 0; tag = 2;
        MPI_Isend(&sd_brasil, 1, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[2]));

        MPI_Wait(&(requests[0]), &(statuss[0])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[1]), &(statuss[1])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[2]), &(statuss[2])); // agora espera terminar a comm coletiva nao bloqueante

        free(mean_cidades);
        free(mean_regioes);
        free(sd_cidades);
        free(sd_regioes);
    }

    if(my_rank == 5){
        MPI_Request requests[3];
        MPI_Status status[3];

        double * median_cidades = (double*) malloc(RCA[0] * RCA[1] * sizeof(double));
        double * median_regioes = (double*) malloc(RCA[0] * sizeof(double));
        double median_brasil;

        #pragma omp parallel for num_threads(T) shared(median_cidades)
        for(int i=0; i<RCA[0] * RCA[1]; i++){
            qsort(&(grades[i * RCA[2]]), RCA[2], sizeof(int), cmpfunc);
            median_cidades[i] = (RCA[2] % 2) ? grades[RCA[2]/2 + i*RCA[2]] : 
                (grades[RCA[2]/2 + i*RCA[2] - 1] + grades[RCA[2]/2 + i*RCA[2]]) / 2.0;

            //printf("Cidade: %d | Median: %.2lf\n",i,median_cidades[i]);
        }

        dst = 0; tag = 0;
        MPI_Isend(median_cidades, RCA[0] * RCA[1], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[0]));

        #pragma omp parallel for num_threads(T) shared(median_regioes)
        for(int i=0; i<RCA[0]; i++){
            qsort(&(grades[i * RCA[1] * RCA[2]]), RCA[1] * RCA[2], sizeof(int), cmpfunc);
            median_regioes[i] = ((RCA[1] * RCA[2]) % 2) ? grades[(i*RCA[1]*RCA[2] + RCA[1]*RCA[2]/2)] : 
                (grades[(i*RCA[1]*RCA[2] + RCA[1]*RCA[2]/2) - 1] + grades[(i*RCA[1]*RCA[2] + RCA[1]*RCA[2]/2)]) / 2.0;
            //printf("Região: %d | Median: %.2lf\n",i,median_regioes[i]);
        }

        dst = 0; tag = 1;
        MPI_Isend(median_regioes, RCA[0], MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[1]));

        qsort(grades, RCA[0] * RCA[1] * RCA[2], sizeof(int), cmpfunc);
        median_brasil = ((RCA[0] * RCA[1] * RCA[2]) % 2) ? grades[((RCA[0] * RCA[1] * RCA[2]))/2] : 
            (grades[((RCA[0] * RCA[1] * RCA[2]))/2 - 1] + grades[((RCA[0] * RCA[1] * RCA[2]))/2]) / 2.0;
        //printf("Brasil | Median: %.2lf\n",median_brasil);

        dst = 0; tag = 2;
        MPI_Isend(&median_brasil, 1, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &(requests[2]));

        MPI_Wait(&(requests[0]), &(status[0])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[1]), &(status[1])); // agora espera terminar a comm coletiva nao bloqueante
        MPI_Wait(&(requests[2]), &(status[2])); // agora espera terminar a comm coletiva nao bloqueante

        free(median_cidades);
        free(median_regioes);
        
    }



    MPI_Finalize();
    return 0;
}
