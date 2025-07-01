#!/bin/bash
# ==========================================================================
# Fast-Matrix-Multiplication library
#   https://github.com/jgdumas/Fast-Matrix-Multiplication
#   Accurate fast matrix multiplications via recursion
# ==========================================================================
# sms2maple: generate maple programs from an L,R,P bilinear algorithm
# ==========================================================================
#   Authors:
#   [J-G. Dumas, C. Pernet, A. Sedoglavic;
#    Strassen's algorithm is not optimally accurate;
#    ISSAC 2024, Raleigh, NC USA, pp. 254-263.
#    https://hal.science/hal-04441653]
# ==========================================================================


# ==========================================================================
##########
# L,R,P matrices
#
Lsms=$1
Rsms=$2
Psms=$3
Suff=check
OPTFLAGS=

if [[ $# -ge 4 ]]; then
   Suff=$4
fi
if [[ $# -ge 5 ]]; then
   OPTFLAGS="-O $5"
fi

Lmat=`dirname $Lsms`/`basename $Lsms .sms`
Rmat=`dirname $Lsms`/`basename $Rsms .sms`
Pmat=`dirname $Lsms`/`basename $Psms .sms`
#echo "$Lmat $Rmat $Pmat"

# ==========================================================================
##########
# functions for program generation
#
source functions4sms.sh



##########
# Extract dimensions
#
Ld=(`grep -v '#' ${Lsms} | head -1 | cut -d' ' -f 1,2`)
Rd=(`grep -v '#' ${Rsms} | head -1 | cut -d' ' -f 1,2`)
Pd=(`grep -v '#' ${Psms} | head -1 | cut -d' ' -f 1,2`)

n=`echo "sqrt(${Rd[1]}*${Pd[0]}/${Ld[1]})"|bc`
m=$(( ${Pd[0]} / ${n} ))
k=$(( ${Rd[1]} / ${n} ))

r=`grep -v '#' ${Lsms} | head -1 | cut -d' ' -f 1`


##########
# test for maple program
#
MAPLEPROG=maple
MAPLEHERE=true
if ! command -v maple &> /dev/null
then
    echo "maple could not be found."
    MAPLEHERE=false
else
    if ! echo "1234+5678;" | maple | grep 6912 &> /dev/null
    then
	echo "`command -v maple` is here but not running"
	MAPLEHERE=false
    else
	echo "maple is up and running"
    fi
fi
if [ "${MAPLEHERE}" != true ] ; then
    MAPLEPROG=tee > /dev/null
fi


# ==========================================================================
# Check MM & Check/Generate associated SLPs
#
sms2slp ${Lmat} ${Rmat} ${Pmat} 1 0 0 "${OPTFLAGS}"

# ==========================================================================
# Bilinear algorithm for matrix multiplication only
#
Lslp=${Lmat}.slp
Rslp=${Rmat}.slp
Pslp=${Pmat}.slp

slp2MMmpl ${Lslp} ${Rslp} ${Pslp} ${m} ${k} ${n} ${r} ${Suff}


# ==========================================================================
# ==========================================================================
