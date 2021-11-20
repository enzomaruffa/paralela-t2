#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define STD_TAG 0
#define ROOT 0

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS 2147483647
#define MULTIPLIER 48271
#define DEFAULT 123456789
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

static long seed = DEFAULT;

double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. 
 * ----------------------------------------------------------------
 */
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
  long t;

  t = MULTIPLIER * (seed % Q) - R * (seed / Q);
  if (t > 0)
    seed = t;
  else
    seed = t + MODULUS;
  return ((double)seed / MODULUS);
}

/*
 * End of the pRNG algorithm
 */

typedef struct
{
  double x, y, z;
  double mass;
} Particle;
typedef struct
{
  double xold, yold, zold;
  double fx, fy, fz;
} ParticleV;

void InitParticles(Particle[], ParticleV[], int);
double ComputeForces(Particle[], Particle[], ParticleV[], int, int, int);
double ComputeNewPos(Particle[], ParticleV[], int, double, int);

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  int rank, proc_count;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc_count);

  MPI_Status status;

  int npart;
  int cnt;      /* number of times in loop */
  double sim_t = 0.0;

  if (rank == 0) {
    int tmp;
    tmp = fscanf(stdin, "%d\n", &npart);
    tmp = fscanf(stdin, "%d\n", &cnt);
  }

  MPI_Bcast(&npart, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
  MPI_Bcast(&cnt, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
  
  // printf("[debug] Process [%d] has the npart and cnt!\n", rank);

  // ========== Defines the indexes

  // MPI_Recv(&npart, 1, MPI_INT, MPI_ANY_SOURCE, STD_TAG, MPI_COMM_WORLD, &status);
  float percentage = (float)rank / proc_count;
  int starting_i = (int)ceil(percentage * npart);
  float next_percentage = (float)(rank+1) / proc_count;
  int end_i = MIN((int)ceil(next_percentage * npart), npart);
  int my_npart = end_i - starting_i;

  // printf("%d: starting_i %d, end_i %d, my_npart %d\n", rank, starting_i, end_i, my_npart);

  int starting_is[proc_count];
  int nparts[proc_count];

  for (int p = 0; p < proc_count; p++) {
    float percentage = (float)p / proc_count;
    int starting_i = (int)ceil(percentage * npart);
    float next_percentage = (float)(p+1) / proc_count;
    int end_i = MIN((int)ceil(next_percentage * npart), npart);

    starting_is[p] = starting_i;
    nparts[p] = end_i - starting_i;
    
    // printf("[debug] [%d] %d: starting_i %d, nparts %d\n", rank, p, starting_is[p], nparts[p]);
  }

  // printf("[debug] Process [%d] has the starting_is and nparts!\n", rank);

  // =========== Instances stuff
  
  Particle *particles; /* Particles */
  ParticleV *pv; /* Particle velocity */
  ParticleV *all_pv; /* Particle velocity for ranks to use as buffer as all */
  
  /* Allocate memory for particles */

  particles = (Particle *)malloc(sizeof(Particle) * npart);
  pv = (ParticleV *)malloc(sizeof(ParticleV) * my_npart);
  all_pv = (ParticleV *)malloc(sizeof(ParticleV) * npart); 

  // =========== Data types
 
  MPI_Datatype MPI_PARTICLE;
  MPI_Type_contiguous(4, MPI_DOUBLE, &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);
  
  MPI_Datatype MPI_PARTICLEV;
  MPI_Type_contiguous(6, MPI_DOUBLE, &MPI_PARTICLEV);
  MPI_Type_commit(&MPI_PARTICLEV);

  if (rank == 0) {
    /* Generate the initial values */
    InitParticles(particles, all_pv, npart);
  }
  
  // rank 0 broadcasts positions of all particles velocities
  MPI_Scatterv(all_pv, nparts, starting_is,
                MPI_PARTICLEV, pv, nparts[rank],
                MPI_PARTICLEV, ROOT, MPI_COMM_WORLD);

  // printf("[debug] Scatter done! Process [%d] now should have the particles velocities\n", rank);
  // if (rank == 0) {
  //   printf("[debug] All pvs: \n");
  //   for (int p = 0; p < npart; p++)
  //     printf("[debug] Rank: %d, [%d]=%f\n", rank, p, all_pv[p].fx);
  // }
  // printf("[debug] PVs: \n");
  // for (int p = 0; p < nparts[rank]; p++)
  //   printf("[debug] Process [%d]: [%d]=%f\n", rank, p, pv[p].fx);

  // rank 0 broadcasts positions of all particles
  MPI_Bcast(particles, npart, MPI_PARTICLE, ROOT, MPI_COMM_WORLD); 
  // printf("[debug] Broadcast done done! Process [%d] now should have the particles\n", rank);
  // printf("[debug] [%d] Particles: \n", rank);
  // for (int p = 0; p < npart; p++)
  //   printf("[debug] Rank: %d, [%d.x]=%f, [%d.y]=%f\n", rank, p, particles[p].x, p, particles[p].y);

  while (cnt--)
  {
    double max_f;
    double final_max_f;
    /* Compute forces (2D only) */

    // each rank computes the forces for their group of particles
    max_f = ComputeForces(particles, particles, pv, npart, starting_i, my_npart);
    // printf("[debug] [%d] PVs: \n", rank);
    // for (int p = 0; p < my_npart; p++)
    //   printf("[debug] Rank: %d, [%d.fx]=%f, [%d.fy]=%f\n", rank, p, pv[p].fx, p, pv[p].fy);

    // MPI_Barrier(MPI_COMM_WORLD);

    // Need to reduce the max_f to all processes
    MPI_Allreduce(&max_f, &final_max_f, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    // printf("[debug] MPI_Allreduce done done! Process [%d] now should have the proper final_max_f=%f\n", rank, final_max_f);

    // other ranks updates the positions 
    /* Once we have the forces, we compute the changes in position */
    sim_t += ComputeNewPos(particles, pv, my_npart, final_max_f, starting_i);
    // printf("[debug] New positions were calculated!\n");
    // printf("[debug] [%d] Particles: \n", rank);
    // for (int p = 0; p < npart; p++)
    //   printf("[debug] Rank: %d, [%d.x]=%f, [%d.y]=%f\n", rank, p, particles[p].x, p, particles[p].y);

    // each rank gathers the particles
    // printf("[debug] [%d] For MPI_Allgatherv, starting address is %p, sending %d\n", rank, particles + nparts[rank], nparts[rank]);

    MPI_Allgatherv(particles + starting_is[rank], nparts[rank], MPI_PARTICLE,
                   particles, nparts, starting_is,
                   MPI_PARTICLE, MPI_COMM_WORLD);

  }

  if (rank == 0) {
    int i = 0;
    for (i = 0; i < npart; i++)
      fprintf(stdout, "%.5lf %.5lf %.5lf\n", particles[i].x, particles[i].y, particles[i].z);
  }
  MPI_Type_free(&MPI_PARTICLE);
  MPI_Type_free(&MPI_PARTICLEV);

  MPI_Finalize();
  return 0;
}

void InitParticles(Particle particles[], ParticleV pv[], int npart)
{
  int i;
  for (i = 0; i < npart; i++)
  {
    particles[i].x = Random();
    particles[i].y = Random();
    particles[i].z = Random();
    particles[i].mass = 1.0;
    pv[i].xold = particles[i].x;
    pv[i].yold = particles[i].y;
    pv[i].zold = particles[i].z;
    pv[i].fx = 0;
    pv[i].fy = 0;
    pv[i].fz = 0;
  }
}

double ComputeForces(Particle myparticles[], Particle others[], ParticleV pv[], int npart, int my_starting_i, int my_npart)
{
  double max_f;
  int i;
  max_f = 0.0;
  for (i = 0; i < my_npart; i++)
  {
    int j;
    double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;

    int my_particle_i = i + my_starting_i;

    rmin = 100.0;
    xi = myparticles[my_particle_i].x;
    yi = myparticles[my_particle_i].y;
    fx = 0.0;
    fy = 0.0;
    for (j = 0; j < npart; j++)
    {
      rx = xi - others[j].x;
      ry = yi - others[j].y;
      mj = others[j].mass;
      r = rx * rx + ry * ry;
      /* ignore overlap and same particle */
      if (r == 0.0)
        continue;
      if (r < rmin)
        rmin = r;
      r = r * sqrt(r);
      fx -= mj * rx / r;
      fy -= mj * ry / r;
    }
    pv[i].fx += fx;
    pv[i].fy += fy;
    fx = sqrt(fx * fx + fy * fy) / rmin;
    if (fx > max_f)
      max_f = fx;
  }
  return max_f;
}

double ComputeNewPos(Particle particles[], ParticleV pv[], int my_npart, double max_f, int my_starting_i)
{
  int i;
  double a0, a1, a2;
  static double dt_old = 0.001, dt = 0.001;
  double dt_new;
  a0 = 2.0 / (dt * (dt + dt_old));
  a2 = 2.0 / (dt_old * (dt + dt_old));
  a1 = -(a0 + a2);
  
  for (int i = 0; i < my_npart; i++)
  {
    int my_particle_i = i + my_starting_i;
    double xi, yi;
    xi = particles[my_particle_i].x;
    yi = particles[my_particle_i].y;
    particles[my_particle_i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
    particles[my_particle_i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
    pv[i].xold = xi;
    pv[i].yold = yi;
    pv[i].fx = 0;
    pv[i].fy = 0;
  }

  dt_new = 1.0 / sqrt(max_f);
  /* Set a minimum: */
  if (dt_new < 1.0e-6)
    dt_new = 1.0e-6;
  /* Modify time step */
  if (dt_new < dt)
  {
    dt_old = dt;
    dt = dt_new;
  }
  else if (dt_new > 4.0 * dt)
  {
    dt_old = dt;
    dt *= 2.0;
  }
  return dt_old;
}
