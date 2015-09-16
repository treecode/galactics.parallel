#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <sys/types.h> /* pid_t */
#include <unistd.h>  /* _exit, fork */
#include <sys/wait.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

extern "C"
{
  void gen_disk(int nptcl, int myseed, float *buffer, int verbose);
  void gen_bulge(int nptcl, int myseed, float *buffer, int verbose);
  void gen_halo (int nptcl, int myseed, float *buffer, int verbose);
}

struct Particle
{
  int id;
  float mass;
  float x,y,z;
  float vx,vy,vz;
};

int childId = 0;
int main(int argc , char * argv[])
{
  const int ndisk  = 1000000;
  const int nbulge =  200000;
  const int nhalo  =  200000;
  
  const int NCHILD = 16;
  const int nptcl = ndisk+nbulge+nhalo;
  const int nByteMap = nptcl*NCHILD*sizeof(Particle);

  /* allocate shared memory */
  pid_t pid;
  int shmfd = shm_open("/SHAREDMEM", O_CREAT|O_RDWR, S_IRUSR|S_IWUSR);
  if (shmfd == -1)
  {
    perror("creator:shm_open");
    exit(EXIT_FAILURE);
  }
  if (ftruncate(shmfd,nByteMap))
  {
    perror("craetor:ftruncate");
    exit(EXIT_FAILURE);
  }

  void *data_ptr = mmap(NULL,nByteMap,PROT_READ|PROT_WRITE, MAP_SHARED,shmfd,0);
  assert(data_ptr != MAP_FAILED);

  for (int i = 0; i < NCHILD; i++)
  {
    pid = fork();
    if (pid < 0)
    {
      fprintf(stderr, " parent :: failed to fork ... \n");
      exit(1);
    } 
    else if (pid == 0)
      break;
    childId++;
  }

  if (pid > 0)
    fprintf(stderr," parent :: done forking ... \n");

  assert(pid >= 0);

#if 0
  if (pid == 0)  
    fprintf(stderr,"I am a child process: id= %d  pid = %d\n", childId, getpid());
  else if (pid  > 0)
    fprintf(stderr,"I am a parent process: pid = %d\n", getpid());
  else
    assert(0);
#endif

  if (pid == 0)
  {
    int shmfd = shm_open("/SHAREDMEM",O_RDWR,0);
    if(shmfd == -1)
    {
      perror("child:shm_open");
      exit(EXIT_FAILURE);
    }
    void *data_ptr = mmap(NULL,nByteMap,PROT_READ|PROT_WRITE, MAP_SHARED,shmfd,0);
    Particle *data = &((Particle*)data_ptr)[nptcl*childId];

    /* generate galaxy */

    const int verbose = childId == 0;
    gen_disk (ndisk,  12345+childId*34, (float*)&data[0           ], verbose);
    gen_bulge(nbulge, 56744+childId*13, (float*)&data[ndisk       ], verbose);
    gen_halo (nhalo,  54837+childId*65, (float*)&data[ndisk+nbulge], verbose);

    return 0;
  }

  if (pid > 0)
  {
    int status;
    for (int i = 0; i < NCHILD; i++)
    {
      int wpid = wait(&status);
//      fprintf(stderr,"Child pid= %d done with status= %d \n", wpid, status);
    }
    fprintf(stderr,"\n parent :: Galaxy generation compleete \n");
    Particle *ptcl = (Particle*)data_ptr;

    /* collect the results from child processes */

    const int ntot = nptcl*NCHILD;
    fprintf(stderr, "nptcl_per_proc= %d  nproc= %d  ntot= %d\n",
        nptcl, NCHILD, ntot);
    const float mscale = 1.0/NCHILD;
#if 0
    for (int i = 0; i < ntot; i++)
    {
      printf(" %.9d  %e  %e %e %e  %e %e %e \n",
          ptcl[i].id,
          ptcl[i].mass*mscale,
          ptcl[i].x,
          ptcl[i].y,
          ptcl[i].z,
          ptcl[i].vx,
          ptcl[i].vy,
          ptcl[i].vz);
    }
#endif
  }

  close(shmfd);
  munmap((void*)data_ptr,nByteMap);
  shm_unlink("/SHAREDMEM");
  return 0;
}
