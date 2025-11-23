@echo off
REM =====================================================================
REM == GROMACS Production Run Script with a FOR Loop for Windows       ==
REM == This script automates the production run in multiple chunks.    ==
REM =====================================================================

REM --- Configuration ---
SET mdp_file=step7_production.mdp
SET last_equi_step=step6.6_equilibration
SET total_chunks=10

REM #####################################################################
REM # CRITICAL: Enable Delayed Expansion to read variables inside a loop #
REM # This allows us to use !variable! to get the current value.         #
REM #####################################################################
SETLOCAL EnableDelayedExpansion

echo Starting GROMACS Production Run Workflow in a loop...
echo.

REM --- Main Loop for Production Chunks ---
REM This loop will run from 1 to %total_chunks%
FOR /L %%i IN (6, 1, %total_chunks%) DO (
    echo [INFO] Preparing Production Chunk %%i of %total_chunks%...

    IF %%i == 1 (
        REM --- This is the FIRST chunk ---
        REM It starts from the last equilibration step and has no checkpoint file (-t)
        echo      -> Special case: Starting from equilibration output.
        gmx grompp -f %mdp_file% -o step7_%%i.tpr -c %last_equi_step%.gro -p topol.top -n index.ndx

    ) ELSE (
        REM --- This is a SUBSEQUENT chunk (2, 3, 4, etc.) ---
        REM It continues from the PREVIOUS chunk using its .gro and .cpt files.
        
        REM Calculate the previous chunk number
        SET /A prev_chunk=%%i-1
        
        echo      -> General case: Continuing from chunk !prev_chunk!.
        gmx grompp -f %mdp_file% -o step7_%%i.tpr -c step7_!prev_chunk!.gro -t step7_!prev_chunk!.cpt -p topol.top -n index.ndx
    )
    
    REM --- Run the simulation for the current chunk ---
    echo      -> Running mdrun for step7_%%i...
    gmx mdrun -v -deffnm step7_%%i
    echo [SUCCESS] Chunk %%i finished.
    echo.
)

echo =================================================================
echo All %total_chunks% production chunks are finished.
echo =================================================================
echo.

pause