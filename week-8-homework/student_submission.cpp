//
// Created by Vincent Bode on 08/07/2020.
//

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "life.h"
#include "Utility.h"
#include "VideoOutput.h"

#include <mpi.h>
#include <omp.h>
#include <iostream>

#define CODE_VERSION 1
/*
  Apply the game of life rules on a Torus --> grid contains shadow rows and columns
  to simplify application of rules i.e. grid actually ranges from grid [ 1.. height - 2 ][ 1 .. width - 2]
*/
void evolve(ProblemData &problemData)
{
    auto &grid = *problemData.readGrid;
    auto &writeGrid = *problemData.writeGrid;
    // For each cell
    for (int i = 1; i < GRID_SIZE - 1; i++)
    {
        for (int j = 1; j < GRID_SIZE - 1; j++)
        {
            // Calculate the number of neighbors
            int sum = grid[i - 1][j - 1] + grid[i - 1][j] + grid[i - 1][j + 1] +
                      grid[i][j - 1] + grid[i][j + 1] +
                      grid[i + 1][j - 1] + grid[i + 1][j] + grid[i + 1][j + 1];

            if (!grid[i][j])
            {
                // If a cell is dead, it can start living by reproduction or stay dead
                if (sum == 3)
                {
                    // reproduction
                    writeGrid[i][j] = true;
                }
                else
                {
                    writeGrid[i][j] = false;
                }
            }
            else
            {
                // If a cell is alive, it can stay alive or die through under/overpopulation
                if (sum == 2 || sum == 3)
                {
                    // stays alive
                    writeGrid[i][j] = true;
                }
                else
                {
                    // dies due to under or overpopulation
                    writeGrid[i][j] = false;
                }
            }
        }
    }
}

/*
  Copies data from the inner part of the grid to
  shadow (padding) rows and columns to transform the grid into a torus.
*/
void copy_edges(bool (&grid)[GRID_SIZE][GRID_SIZE])
{
    // Copy data to the boundaries
    for (int i = 1; i < GRID_SIZE - 1; i++)
    {
        // join rows together
        grid[i][0] = grid[i][GRID_SIZE - 2];
        grid[i][GRID_SIZE - 1] = grid[i][1];
    }

    for (int j = 1; j < GRID_SIZE - 1; j++)
    {
        // join columns together
        grid[0][j] = grid[GRID_SIZE - 2][j];
        grid[GRID_SIZE - 1][j] = grid[1][j];
    }

    // Fix corners
    grid[0][0] = grid[GRID_SIZE - 2][GRID_SIZE - 2];
    grid[GRID_SIZE - 1][GRID_SIZE - 1] = grid[1][1];
    grid[0][GRID_SIZE - 1] = grid[GRID_SIZE - 2][1];
    grid[GRID_SIZE - 1][0] = grid[1][GRID_SIZE - 2];
}

int main(int argc, char** argv) {

    MPI_Init (& argc , & argv ); /* starts MPI */

    int rank , size;

    
    MPI_Comm_rank ( MPI_COMM_WORLD , & rank ); /* process id */
    MPI_Comm_size ( MPI_COMM_WORLD , & size ); /* number processes */
    
    // As with Jack Sparrow's exercise, this needs FFMPEG (new and improved: this now works with more video players).
    // As an alternative, you can write individual png files to take a look at the data.
    //if (rank == 0) {
        bool activateVideoOutput = false;
        if (argc > 1) {
            if (argc == 2 && strcmp(argv[1], "-g") == 0) {
                activateVideoOutput = true;
            }
            else {
                fprintf(stderr, "Usage:\n  %s [-g]\n    -g: Activate graphical output.\n", argv[0]);
                exit(153);
            }
        }

        auto* problemData = new ProblemData;

        if (activateVideoOutput) {
            VideoOutput::beginVideoOutput();
            VideoOutput::saveToFile(*problemData->readGrid, "grid_beginning.png");
        }

        
        if (rank == 0) {
            Utility::readProblemFromInput(CODE_VERSION, *problemData);
        }

    // the number of rows calculated by a process
    int proc_rows = GRID_SIZE / size;

    // the starting row in this process
    int st_row = rank * proc_rows;

    // the end row(not included) in this process
    int end_row = (rank + 1) * proc_rows;
    //std::cout << end_row << "   " <<size << std::endl;
   
    //bool *temp = &(*problemData->readGrid)
    MPI_Bcast(problemData->readGrid, (GRID_SIZE) * (GRID_SIZE), MPI_C_BOOL, 0, MPI_COMM_WORLD);
    //for (int i = 0; i < GRID_SIZE; ++i){
     //   MPI_Bcast(problemData->readGrid, (GRID_SIZE), MPI_FLOAT, 0, MPI_COMM_WORLD);
    //}
    //std::cout << "DONE2" << std::endl;
    
    //#pragma omp parallel for schedule(dynamic)
    //TODO@Students: This is the main simulation. Parallelize it using MPI.
    for (int iteration = 0; iteration < NUM_SIMULATION_STEPS; ++iteration) {
        copy_edges(*problemData->readGrid);
        auto& grid = *problemData->readGrid;
      //for each processor
      for(int i = 0; i < size; ++i){
        for(int j = 0; j < buffer_length; ++j){
            buffer_top[j] = grid[][];
            buffer_bottom[j] = grid[][];
        } 
        if (rank == 0)
        {
            std::cout << "DONE1" << std::endl;
            // send lowest row in block to process with rank 1
            MPI_Send(&grid[end_row - 1][0], GRID_SIZE, MPI_C_BOOL, 1, 0, MPI_COMM_WORLD);
            // MPI_Send((problemData + (end_row - 1) * GRID_SIZE)->readGrid, GRID_SIZE, MPI_C_BOOL, 1, 0, MPI_COMM_WORLD);
           
            //(problemData + GRID_SIZE)->readGrid
            //std::cout << "DONE2" << std::endl;
            // receive upper ghost layer from process with rank size-1
            MPI_Recv(*problemData->readGrid + (GRID_SIZE - 1) * GRID_SIZE, GRID_SIZE, MPI_FLOAT, size - 1, 0, MPI_COMM_WORLD, nullptr);

            // send uppermost row in block to process with rank size-1
            MPI_Send(*problemData->readGrid + st_row * GRID_SIZE, GRID_SIZE, MPI_FLOAT, size - 1, 1, MPI_COMM_WORLD);

            // receive lower ghost layer from process with rank 1
            std::cout << "DONE3" << std::endl;
            MPI_Recv(problemData->readGrid[end_row], GRID_SIZE, MPI_C_BOOL, 1, 1, MPI_COMM_WORLD, nullptr);
            std::cout << "DONE4" << std::endl;

        }
        else if (rank == size - 1)
        {
            // receive upper ghost layer from process with rank size - 2
            MPI_Recv(problemData->readGrid[st_row - 1], GRID_SIZE, MPI_C_BOOL, rank - 2, 0, MPI_COMM_WORLD, nullptr);

            // send uppermost row in block to process with rank size - 2
            MPI_Send(problemData->readGrid[st_row], GRID_SIZE, MPI_C_BOOL, rank - 2, 1, MPI_COMM_WORLD);
        }
        else
        {
            // TODO: receive upper ghost layer from process with rank rank-1
            MPI_Recv(problemData->readGrid[st_row - 1], GRID_SIZE, MPI_C_BOOL, rank-1, 0, MPI_COMM_WORLD, nullptr);
            // TODO: send lowest row in block to the process with rank (rank+1)%size
            MPI_Send(problemData->readGrid[end_row - 1], GRID_SIZE, MPI_C_BOOL, (rank+1)%size, 0, MPI_COMM_WORLD);
            // TODO: receive lower ghost layer from process with rank (rank+1)%size
            MPI_Recv(problemData->readGrid[end_row % GRID_SIZE], GRID_SIZE, MPI_C_BOOL, (rank+1)%size, 1, MPI_COMM_WORLD, nullptr);
            // TODO: send uppermost row in block to the process with rank rank-1
            MPI_Send(problemData->readGrid[st_row], GRID_SIZE, MPI_C_BOOL, rank-1, 1, MPI_COMM_WORLD);
        }
      }


        if (rank == 0) {
            if (activateVideoOutput) {
                VideoOutput::writeVideoFrames(*problemData);
            }

            if (iteration % SOLUTION_REPORT_INTERVAL == 0) {
                Utility::outputIntermediateSolution(iteration, *problemData);
            }
        }
        evolve(*problemData);

        problemData->swapGrids();

    }
    //delete problemData;
    
    MPI_Finalize();
    return 0;
}
