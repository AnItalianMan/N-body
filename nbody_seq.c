#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

#define SOFTENING 1e-9f

typedef struct
{
    float x, y, z, vx, vy, vz;
} Body;

void randomizeBodies(float *data, int n)
{
    int multiplier[] = {2, 1, 4, 1};
    for (int i = 0; i < n; i++)
    {
        data[i] = 2.0f * (i / 50) + 1.0f * multiplier[i % 4];
    }
}

void bodyForce(Body *p, float dt, int n)
{
    for (int i = 0; i < n; i++)
    {
        float Fx = 0.0f;
        float Fy = 0.0f;
        float Fz = 0.0f;

        for (int j = 0; j < n; j++)
        {
            float dx = p[j].x - p[i].x;
            float dy = p[j].y - p[i].y;
            float dz = p[j].z - p[i].z;
            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }

        p[i].vx += dt * Fx;
        p[i].vy += dt * Fy;
        p[i].vz += dt * Fz;
    }
}

void writeResults(Body *bodies, int size)
{
    FILE *fp;

    fp = fopen("result.txt", "w");

    fprintf(fp, "---------------------- Processed %d bodies ----------------------\n", size);
    for (int i = 0; i < size; i++)
    {
        fprintf(fp, "+ body number: [%d]\n", i);
        fprintf(fp, "|--  x: %f     y: %f     z: %f\n", bodies[i].x, bodies[i].y, bodies[i].z);
        fprintf(fp, "|-- vx: %f    vy: %f    vz: %f\n", bodies[i].vx, bodies[i].vy, bodies[i].vz);
        fprintf(fp, "-----------------------------------------------------------------------------------------\n");
    }
    fclose(fp);
}

int main(const int argc, const char **argv)
{

    int nBodies = 30;

    const float dt = 0.01f;
    const int nIters = 25;

    int bytes = nBodies * sizeof(Body);
    float *buf = (float *)malloc(bytes);
    Body *p = (Body *)buf;

    randomizeBodies(buf, 6 * nBodies);

    double totalTime = 0.0;

    for (int iter = 1; iter <= nIters; iter++)
    {
        bodyForce(p, dt, nBodies);

        for (int i = 0; i < nBodies; i++)
        {
            p[i].x += p[i].vx * dt;
            p[i].y += p[i].vy * dt;
            p[i].z += p[i].vz * dt;
        }
    }

    writeResults(p, nBodies);
    free(buf);
}
