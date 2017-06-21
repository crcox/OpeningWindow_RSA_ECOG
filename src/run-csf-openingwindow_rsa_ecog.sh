#!/bin/bash

# -- SGE optoins (whose lines must begin with #$)

#$ -S /bin/bash # The jobscript is written for the bash shell
#$ -V # Inherit environment settings (e.g., from loaded modulefiles)
#$ -o ./logfiles
#$ -e ./logfiles
#$ -cwd # Run the job in the current directory

# THIS MUST BE UPDATED TO REFLECT THE NUMBER OF JOBS! THESE DO NOT NEED TO
# MATCH THE DIRECTORY NAMES!
#$ -t 1-80

# Relevant modules? Perhaps?
# compilers/gcc/4.6.2
# compilers/gcc/4.7.0 * most current that is officially compatible with mex
# compilers/gcc/4.8.2
# compilers/gcc/4.9.0
# compilers/gcc/6.3.0
#
# apps/binapps/matlab/R2010a
# apps/binapps/matlab/R2011a
# apps/binapps/matlab/R2012a
# apps/binapps/matlab/R2013a
# apps/binapps/matlab/R2014a
# apps/binapps/matlab/R2015a
# apps/binapps/matlab/R2015aSP1
# -- the commands to be executed (programs to be run) on a compute node:

module load apps/binapps/matlab/R2015aSP1
export MATLABDIR="$(dirname $(dirname $(which matlab)))"

NJOBS=80
NDIGITS=${#NJOBS}
# Subdirectories will all have this common root (saves me some typing)
BASE="$HOME/scratch/results/Naming_ECOG/OpeningWindow/all"

# Path to my executable
EXE="$BASE/shared/run_OpeningWindow_RSA_ECOG.sh"

# Path too the root of the data directory tree
DATAROOT="/mnt/sw01-home01/mbmhscc4/scratch/data/Naming_ECoG/avg"

# Which subset of examples to consider (all, animate, or inanimate)
SUBSET="all"

# Here comes some fancy bash: apparently when you combine two or more brace
# expansions, you get the cartesian product of the two. That is:
#    {A,B}{1,2}
# Will produce:
#    A1
#    A2
#    B1
#    B2
#
# N.B. The right-most brace will be the inner-most loop during the combined
# expansion (A1, A2, then B1, B2).
#
# If you add a space between the two brace expressions, a space will be be
# included between the outputs. Therefore:
#    {A,B}" "{1,2}
# Will produce:
#    A 1
#    A 2
#    B 1
#    B 2
#
# This is convenient, because the space is used to separate array elements.
#
# In the following, the outer parenthesis are critical. It makes sure the
# evalation happens before assignment, and it encodes the expanded output as an
# array. Each element in the array is a pair of values. To continue the example
# from above:
#
#    x=({A,B}" "{1,2})
# Will permit the following commands and results (which are a bit indented):
#    echo ${x[0]}
#      A 1
#    echo ${x[3]}
#      B 2
#    y=(${x[1]})
#    echo ${y[0]}
#      A
#    echo ${y[1]}
#      2
CONDITIONS=({1,2,3,5,7,8,9,10}" "{1..10})
#           Subjects             Holdouts

echo "Number of jobs: $NJOBS"

# Remember that $SGE_TASK_ID will be 1, 2, 3, ... 24.
# BASH array indexing starts from zero so decrment.
TID=$[SGE_TASK_ID-1]

# Index in to the arrays of directory names to create a path
COND=(${CONDITIONS[$TID]})
SUBJ=${COND[0]}
HOLD=${COND[1]}

# Echo some info to the job output file
echo "Running SGE_TASK_ID $SGE_TASK_ID, Subject ${SUBJ} and Holdout ${HOLD}"

# Finally run my executable from the correct directory
$EXE "$MATLABDIR" "$SUBSET" $SUBJ $HOLD "parallel" 0 "dataroot" "$DATAROOT"
