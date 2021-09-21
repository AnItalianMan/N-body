#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#define SOFTENING 1e-9f
#define MASTER 0

typedef struct
{
    float x, y, z, vx, vy, vz;

} Body;

void randomizeBodies(float *data, int n)
{
    int multiplier[] = {2, 1, 4, 1};
    for (int i = 0; i < n; i++)
    {
        data[i] = 2.0f * (i / (float)50) + 1.0f * multiplier[i % 4];
    }
}

void applyBodyForce(Body *myBody, int myBodySize, Body *otherBodies, int otherBodiesSize, float dt)
{
    for (int i = 0; i < myBodySize; i++)
    {
        float Fx = 0.0f;
        float Fy = 0.0f;
        float Fz = 0.0f;

        for (int j = 0; j < otherBodiesSize; j++)
        {
            float dx = otherBodies[j].x - myBody[i].x;
            float dy = otherBodies[j].y - myBody[i].y;
            float dz = otherBodies[j].z - myBody[i].z;

            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }

        myBody[i].vx += dt * Fx;
        myBody[i].vy += dt * Fy;
        myBody[i].vz += dt * Fz;
    }
}

void initDisplsAndSendCounts(int processes, int nBodies, int *displacements, int *sendCounts, int minSize)
{
    int numberOfProcessorsToUse = nBodies - minSize * processes > 0 ? processes : nBodies / minSize;

    if (numberOfProcessorsToUse == 0)
        numberOfProcessorsToUse = 1;

    int bodiesForEachProcess = nBodies / numberOfProcessorsToUse;
    int resto = nBodies % numberOfProcessorsToUse;

    for (int i = 0; i < numberOfProcessorsToUse; i++)
    {
        sendCounts[i] = bodiesForEachProcess;
        if (resto > 0)
        {
            sendCounts[i]++;
            resto--;
        }
    }

    displacements[0] = 0;
    for (int i = 1; i < processes; i++)
    {
        displacements[i] = displacements[i - 1] + sendCounts[i - 1];
    }
}

Body *getBodyAddressForActualIndex(Body *basicBodyAddress, int myRank, int index, int *sendCounts)
{
    int pos = 0;
    for (int i = 0; i < index; i++)
    {
        if (i != myRank)
            pos += sendCounts[i];
    }
    return basicBodyAddress + pos;
}

void createBodyStruct(MPI_Datatype *bodyType)
{
    MPI_Datatype oldtypes[1];
    int blockcounts[1];
    MPI_Aint offsets[1];

    offsets[0] = 0;
    blockcounts[0] = 6;
    oldtypes[0] = MPI_FLOAT;

    MPI_Type_create_struct(1, blockcounts, offsets, oldtypes, bodyType);
    MPI_Type_commit(bodyType);
}

void initVariables(int nBodies, int *myRank, int *processes, int **sendCounts, int **displacements, Body **recv, int *myNumberOfElements, int minSize)
{
    MPI_Comm_size(MPI_COMM_WORLD, processes);
    MPI_Comm_rank(MPI_COMM_WORLD, myRank);

    *sendCounts = (int *)malloc(*processes * sizeof(int));
    *displacements = (int *)malloc(*processes * sizeof(int));
    memset(*sendCounts, 0, *processes * sizeof(int));

    initDisplsAndSendCounts(*processes, nBodies, *displacements, *sendCounts, minSize);
    *myNumberOfElements = sendCounts[0][*myRank];
    *recv = (Body *)malloc(*myNumberOfElements * sizeof(Body));
}

void getAndSendBodiesToOtherProcesses(Body *otherBodiesBasicAddress, MPI_Datatype bodyType, MPI_Request *requests, int myRank, int processes, int *sendCounts, int myNumberOfElements, Body *recv)
{
    for (int i = 0; i < processes; i++)
    {
        if (i != myRank)
        {
            MPI_Ibcast(otherBodiesBasicAddress, sendCounts[i], bodyType, i, MPI_COMM_WORLD, &requests[i]);
            otherBodiesBasicAddress += sendCounts[i];
        }
        else
        {
            MPI_Ibcast(recv, myNumberOfElements, bodyType, myRank, MPI_COMM_WORLD, &requests[i]);
        }
    }
}

void waitForOtherBodiesAndApplybodyForce(MPI_Request *requests, Body *otherProcessesBodies, Body *recv, int *sendCounts, int myRank, int processes, int myNumberOfElements, float dt)
{
    int i = 0;
    while (i != processes)
    {
        int index;
        MPI_Status status;
        MPI_Waitany(processes, requests, &index, &status);
        i++;
        if (index != myRank)
        {
            Body *actualBody = getBodyAddressForActualIndex(otherProcessesBodies, myRank, index, sendCounts);
            applyBodyForce(recv, myNumberOfElements, actualBody, sendCounts[index], dt);
        }
    }
}

void integrateBodyPosition(Body *recv, int myNumberOfElements, float dt)
{
    for (int i = 0; i < myNumberOfElements; i++)
    {
        recv[i].x += recv[i].vx * dt;
        recv[i].y += recv[i].vy * dt;
        recv[i].z += recv[i].vz * dt;
    }
}

void runSingleIteration(Body *otherProcessesBody, Body *recv, MPI_Datatype bodyType, MPI_Request *requests, int *sendCounts, int myRank, int processes, int myNumberOfElements, float dt)
{
    getAndSendBodiesToOtherProcesses(otherProcessesBody, bodyType, requests, myRank, processes, sendCounts, myNumberOfElements, recv);
    applyBodyForce(recv, myNumberOfElements, recv, myNumberOfElements, dt);
    waitForOtherBodiesAndApplybodyForce(requests, otherProcessesBody, recv, sendCounts, myRank, processes, myNumberOfElements, dt);
    integrateBodyPosition(recv, myNumberOfElements, dt);
}

void runAllIterations(Body *recv, int *sendCounts, MPI_Datatype bodyType, int nBodies, int myNumberOfElements, int myRank, int processes, int nIters)
{
    const float dt = 0.01f;
    int currentIteration = 0;
    while (currentIteration < nIters)
    {
        Body *otherProcessesBody = (Body *)malloc((nBodies - myNumberOfElements) * sizeof(Body));
        MPI_Request *requests = (MPI_Request *)malloc(processes * sizeof(MPI_Request));

        runSingleIteration(otherProcessesBody, recv, bodyType, requests, sendCounts, myRank, processes, myNumberOfElements, dt);

        free(otherProcessesBody);
        free(requests);
        currentIteration++;
    }
}

void initProcesses(Body *recv, int *sendCounts, int *displacements, MPI_Datatype bodyType, int nBodies, int myRank, int myNumberOfElements)
{
    Body *p = NULL;
    if (myRank == MASTER)
    {
        int bytes = nBodies * sizeof(Body);
        float *buf = (float *)malloc(bytes);
        p = (Body *)buf;
        randomizeBodies(buf, 6 * nBodies);
    }

    MPI_Scatterv(p, sendCounts, displacements, bodyType,
                 recv, myNumberOfElements, bodyType, MASTER, MPI_COMM_WORLD);
    free(p);
}

void writeResults(Body *bodies, int size, int myRank, double time)
{
    if (myRank == MASTER)
    {
        FILE *fp;

        fp = fopen("result.txt", "w");

        fprintf(fp, "---------------------- Processed %d bodies in %f seconds ----------------------\n", size, time);
        for (int i = 0; i < size; i++)
        {
            fprintf(fp, "+ body number: [%d]\n", i);
            fprintf(fp, "|--  x: %f     y: %f     z: %f\n", bodies[i].x, bodies[i].y, bodies[i].z);
            fprintf(fp, "|-- vx: %f    vy: %f    vz: %f\n", bodies[i].vx, bodies[i].vy, bodies[i].vz);
            fprintf(fp, "-----------------------------------------------------------------------------------------\n");
        }
        fclose(fp);
        printf("Execution time: %f\n", time);
    }
}

void collectIterationsResult(Body **processedBodies, Body *recv, int *sendCounts, int *displacements, MPI_Datatype bodyType, int nBodies, int myRank, int myNumberOfElements)
{
    *processedBodies = NULL;

    if (myRank == MASTER)
    {
        *processedBodies = (Body *)malloc(nBodies * sizeof(Body));
    }

    MPI_Gatherv(recv, myNumberOfElements, bodyType, *processedBodies, sendCounts, displacements, bodyType, MASTER, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    const int nBodies = argv[1] != NULL ? atoi(argv[1]) : 30000;
    const int nIters = argv[2] != NULL ? atoi(argv[2]) : 10;
    const int minSize = argv[3] != NULL ? atoi(argv[3]) : 100;

    MPI_Datatype bodyType;
    int myRank;
    int processes;
    int *sendCounts;
    int *displacements;
    int myNumberOfElements;
    Body *recv;
    Body *processedBodies;
    double start, end;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    createBodyStruct(&bodyType);
    initVariables(nBodies, &myRank, &processes, &sendCounts, &displacements, &recv, &myNumberOfElements, minSize);
    initProcesses(recv, sendCounts, displacements, bodyType, nBodies, myRank, myNumberOfElements);

    runAllIterations(recv, sendCounts, bodyType, nBodies, myNumberOfElements, myRank, processes, nIters);

    collectIterationsResult(&processedBodies, recv, sendCounts, displacements, bodyType, nBodies, myRank, myNumberOfElements);
    MPI_Type_free(&bodyType);
    free(recv);
    free(sendCounts);
    free(displacements);

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    writeResults(processedBodies, nBodies, myRank, end - start);
    free(processedBodies);

    MPI_Finalize();
    return 0;
}
