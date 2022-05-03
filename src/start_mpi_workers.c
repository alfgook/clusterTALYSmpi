/* manager */ 
#include "mpi.h"
#include "unistd.h"
#include <R.h>

#define MIN(x,y) ((x<y)?x:y)

void initalize_mpi(int *number_of_workers) {
   int world_size,flag;

   MPI_Init(NULL, NULL);
   printf("mpi session initalized\n");
   printf("number_of_workers = %d\n",*number_of_workers);
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_UNIVERSE_SIZE, &number_of_workers, &flag);
   *number_of_workers = *number_of_workers;
   printf("number_of_workers = %d\n",*number_of_workers);
}

void finalize_mpi() {
   printf("finalize_mpi is called...\n");
   int init_flag;
   MPI_Initialized(&init_flag);
   if(init_flag) {
      MPI_Finalize();
      printf("the mpi session has been finalized\n");
   }
}

void start_mpi_workers(const char **worker_program ,
                        char *job_list[],
                        const int *number_of_jobs,
                        const int *number_of_workers,
                        char **talys_exe,
                        char **bin_path) 
{ 

   //MPI_Init(NULL, NULL);  

   int world_size, universe_size, *universe_sizep, flag; 
   MPI_Comm everyone;           /* intercommunicator */ 

   int nbr_of_workers = (*number_of_workers);
   int nbr_of_jobs = (*number_of_jobs);

   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   if (world_size != 1) {
      printf("Top heavy with management");
      return;
   }

   //MPI_Attr_get(MPI_COMM_WORLD, MPI_UNIVERSE_SIZE, &universe_sizep, &flag);
   MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_UNIVERSE_SIZE, &universe_sizep, &flag);

   if (!flag) { 
      printf("This MPI does not support UNIVERSE_SIZE."); 
      return;
   } else {
      universe_size = *universe_sizep;
   }

   //printf("universe_size = %d\n",universe_size);
   if (universe_size == 1) {
      printf("No room to start workers");
      return;
   }

   if(nbr_of_workers>0) {
      nbr_of_workers = MIN(nbr_of_workers,universe_size-1);
      nbr_of_workers = MIN(nbr_of_workers,nbr_of_jobs);
   } else {
      nbr_of_workers = MIN(universe_size-1,nbr_of_jobs);
   }
   
   //printf("nbr_of_workers = %d\n",nbr_of_workers);
   //printf("nbr_of_jobs = %d\n",nbr_of_jobs);
   //printf("universe_size = %d\n",universe_size);
   //printf("job_list[0] = %s\n",job_list[0]);

   // copy the list of jobs from the calling R-script and 
   char **argv = (char **) malloc((nbr_of_jobs+2) * sizeof(char *));
   for(int i=0;i<nbr_of_jobs;i++) {
      argv[i] = job_list[i];
   }
   argv[nbr_of_jobs] = talys_exe[0];
   if(talys_exe[0]==NULL) argv[nbr_of_jobs] = "talys";
   // add a NULL pointer at the end
   // this is required by MPI_Comm_spawn to calculate the number of arguments
   argv[nbr_of_jobs+1] = NULL;


   // create info object
   MPI_Info info;
   MPI_Info_create(&info);
   // add /usr/local/bin to the info object
   MPI_Info_set(info,"wdir",bin_path[0]);
   // apparently I can only have one value for wdir, so it needs to be known at compile time

   MPI_Comm_spawn(worker_program[0],
                     argv,
                     nbr_of_workers,
                     info, 
                     0,
                     MPI_COMM_SELF,
                     &everyone,
                     MPI_ERRCODES_IGNORE);


   /*MPI_Comm_spawn(worker_program[0],
                     argv,
                     nbr_of_workers,
                     MPI_INFO_NULL, 
                     0,
                     MPI_COMM_SELF,
                     &everyone,
                     MPI_ERRCODES_IGNORE);*/
   
   // Could add parallel code here. The communicator "everyone" can be used to communicate with
   // the spawned processes, which have ranks 0,.. nbr_of_workers-1 in the remote group of 
   // the intercommunicator "everyone".
   // I need to add a wait for the workers to finish thier tasks here!!!
   //MPI_Barrier(everyone);

   printf("%d workers spawned... waiting...\n",nbr_of_workers);
   for(int worker=0;worker<nbr_of_workers;worker++) {
      int rank;
      MPI_Recv(&rank, 1, MPI_INT, worker, 0, everyone, MPI_STATUS_IGNORE);
      //printf("master received %d from worker %d\n",rank,worker);
   }
   

   free(argv);

   printf("all done!\n");
   //MPI_Finalize();
} 