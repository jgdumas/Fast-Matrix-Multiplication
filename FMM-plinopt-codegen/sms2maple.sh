#!/bin/bash
# ==========================================================================
# Fast-Matrix-Multiplication library
#   https://github.com/jgdumas/Fast-Matrix-Multiplication
#   Accurate fast matrix multiplications via recursion
# ==========================================================================
# sms2maple: generate maple programs from an L,R,P bilinear algorithm
# exeample: `echo n | sms2maple.sh data/4x4x4_48_rational_{L,R,P}.sms check`
# ==========================================================================
#   Authors:
#   [J-G. Dumas, C. Pernet, A. Sedoglavic;
#    Strassen's algorithm is not optimally accurate;
#    ISSAC 2024, Raleigh, NC USA, pp. 254-263.
#    https://hal.science/hal-04441653]
# ==========================================================================

##########
# Parsing args
#
MATS=()
MMCHECK=1
Suff=check
OPTFLAGS=""
SQRT=0
REPL=0
EXPO=0
MODU=0
MODP=0
IRRED="Empty"

while [[ $# -gt 0 ]]; do
  case $1 in
    -O|--Optflags)
      OPTFLAGS="-O $2"
      shift # past argument
      shift # past value
      ;;
    -m|-mm|-mmc|--MMcheck)
      MMCHECK=1
      shift # past argument
      ;;
    -n|-nm|-nmm|-nmmc|-nc|--NoMMcheck)
      MMCHECK=0
      shift # past argument
      ;;
    -p|-pm|--PMcheck)
      MMCHECK=0
      PMCHECK=1
      shift # past argument
      ;;
    -I|--Ireducible)
      MODP=1
      IRRED="$2"
      MMCHECK=0
      PMCHECK=1
      shift # past argument
      shift # past value
      ;;
    -q|--prime)
      MODU="$2"
      shift # past argument
      shift # past value
      ;;
    -r|--modular)
      REPL="$2"
      EXPO="$3"
      SQRT="$4"
      shift # past argument
      shift # past value
      shift # past value
      shift # past value
      ;;
    -h|--h|-help|--help|-*|--*)
      echo "Usage: $0 [-O #|-r # # #|-m|-n] L.sms R.sms P.sms [name]"
      echo "  generates a Maple program from L,R,P matrices."
      echo "  [name]: the program suffix (default '${Suff}')."
      echo "  -m/-n: L,R,P are checked/not checked as a mat. mul. (default no)."
      echo "  -O N: optimizer with N loops (default is ${OPTFLAGS})."
      echo "  -r r e s: compute modulo (r^e-s) in sms files (default none)."
      exit 1
      ;;
    *)
      MATS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done
set -- "${MATS[@]}" # restore positional parameters

# ==========================================================================
##########
# L,R,P matrices
#
Lsms=$1
Rsms=$2
Psms=$3
if [[ $# -ge 4 ]]; then
   Suff=$4
fi

Lmat=`dirname $Lsms`/`basename $Lsms .sms`
Rmat=`dirname $Lsms`/`basename $Rsms .sms`
Pmat=`dirname $Lsms`/`basename $Psms .sms`
#echo "$Lmat $Rmat $Pmat"

# ==========================================================================
##########
# functions for program generation
#
source `dirname $0`/functions4sms.sh

##########
# test for PLinOpt/FMM program
#
PLinOptPresent
FMMcodegenPresent
FMMD=`dirname $0`


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
	echo "# Maple is up and running"
    fi
fi
if [ "${MAPLEHERE}" != true ] ; then
    MAPLEPROG=tee > /dev/null
fi


# ==========================================================================
# Check MM & Check/Generate associated SLPs
#
sms2slp ${Lmat} ${Rmat} ${Pmat} ${MMCHECK} ${REPL} ${EXPO} ${SQRT} "${OPTFLAGS}" ${MODU}

# ==========================================================================
# Bilinear algorithm for matrix multiplication only
#
Lslp=${Lmat}.slp
Rslp=${Rmat}.slp
Pslp=${Pmat}.slp


##########
# Extract dimensions
#
Ld=(`grep -v '#' ${Lsms} | head -1 | cut -d' ' -f 1,2`)
Rd=(`grep -v '#' ${Rsms} | head -1 | cut -d' ' -f 1,2`)
Pd=(`grep -v '#' ${Psms} | head -1 | cut -d' ' -f 1,2`)

if [[ "$PMCHECK" -eq 1 ]]; then
    m=${Ld[1]}
    k=${Rd[1]}
    n=${Pd[0]}
    r=${Ld[0]}

    slp2PMmpl ${Lslp} ${Rslp} ${Pslp} ${m} ${k} ${n} ${r} ${Suff} ${MODP} "${IRRED}"
else
    n=`echo "sqrt(${Rd[1]}*${Pd[0]}/${Ld[1]})"|bc`
    m=$(( ${Pd[0]} / ${n} ))
    k=$(( ${Rd[1]} / ${n} ))
    r=`grep -v '#' ${Lsms} | head -1 | cut -d' ' -f 1`

    slp2MMmpl ${Lslp} ${Rslp} ${Pslp} ${m} ${k} ${n} ${r} ${Suff} ${REPL} ${EXPO} ${SQRT}
fi



# ==========================================================================
# ==========================================================================
