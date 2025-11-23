@echo off
REM =================================================================
REM == GROMACS Simulation Script for Windows (Batch File)          ==
REM == This script automates the Minimization and Equilibration    ==
REM =================================================================

echo Starting GROMACS Simulation Workflow...
echo.

REM --- Step 6.6: Equilibration ---
echo [INFO] Step 6.6: Running Equilibration...
gmx grompp -f step6.6_equilibration.mdp -c step6.5_equilibration.gro -r step5_input.gro -p topol.top -n index.ndx -o step6.6_equilibration.tpr
gmx mdrun -v -deffnm step6.6_equilibration
echo [SUCCESS] Step 6.6 finished. Equilibration phase is complete!
echo.

echo ============================================================
echo All equilibration steps are finished.
echo You can now proceed to the production run (step 7).
echo ============================================================
echo.

pauses