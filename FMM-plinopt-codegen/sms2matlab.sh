#!/bin/bash
# ==========================================================================
# Fast-Matrix-Multiplication library
#   https://github.com/jgdumas/Fast-Matrix-Multiplication
#   Accurate fast matrix multiplications via recursion
# ==========================================================================
#   Authors:
#   [J-G. Dumas, C. Pernet, A. Sedoglavic;
#    Strassen's algorithm is not optimally accurate;
#    ISSAC 2024, Raleigh, NC USA, pp. 254-263.
#    https://hal.science/hal-04441653]
# ==========================================================================


##########
# Parsing args
MATS=()
OPTFLAGS="-E -N"
MMCHECK=1
CoBTYPE=0
SQRT=0
PLACE=0
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
    -c|-cob|--CoB)
      MMCHECK=0
      CoBTYPE=1
      shift # past argument
      ;;
    -n|-nm|-nmm|-nmmc|-nc|--NoMMcheck)
      MMCHECK=0
      shift # past argument
      ;;
    -p|--placeholder)
      SQRT="$2"
      PLACE="$3"
      shift # past argument
      shift # past value
      shift # past value
      ;;
    -h|--h|-help|--help|-*|--*)
      echo "Usage: sms2matlab.sh [-O #|-p # #|-m|-n|-c] L.sms R.sms P.sms name"
      echo "  generates matlab program name.m from L,R,P matrices."
      echo "  -m: L,R,P are a matrix multiplication algorithm (default)."
      echo "  -c: L,R,P are change of bases matrices (default MM)."
      echo "  -m/-n: L,R,P are checked/not checked as a mat. mul. (default no)."
      echo "  -O N: optimizer with N loops (default is ${OPTFLAGS})."
      echo "  -p S P: P replaces sqrt(S) in sms files (default none)."
      exit 1
      ;;
    *)
      MATS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done
set -- "${MATS[@]}" # restore positional parameters

##########
# L,R,P matrices
Lsms=$1
Rsms=$2
Psms=$3
File=$4

##########
# Extract dimensions

Ld=(`grep -v '#' ${Lsms} | head -1 | cut -d' ' -f 1,2`)
Rd=(`grep -v '#' ${Rsms} | head -1 | cut -d' ' -f 1,2`)
Pd=(`grep -v '#' ${Psms} | head -1 | cut -d' ' -f 1,2`)

n=`echo "sqrt(${Rd[1]}*${Pd[0]}/${Ld[1]})"|bc`
m=$(( ${Pd[0]} / ${n} ))
k=$(( ${Rd[1]} / ${n} ))

# m=`echo $mkn | sed 's/x/ /g' | cut -d' ' -f4`
# k=`echo $mkn | sed 's/x/ /g' | cut -d' ' -f5`
# n=`echo $mkn | sed 's/x/ /g' | cut -d' ' -f6`
r=`grep -v '#' ${Lsms} | head -1 | cut -d' ' -f 1`


##########
# Generator of SLPs from 3 (2 direct and 1 transposed) sms matrices
function sms2slp {
    local Lsms=$1
    local Rsms=$2
    local Psms=$3
    local MMCHECK=$4
    local SQRT=$5
    local PLACE=$6

##########
# Do represent a matrix multiplication algorithm

if [[ "$MMCHECK" -eq 1 ]]; then
  local MMFLAGS=""
  if [[ "$SQRT" -ne 0 ]]; then
    local MMFLAGS="128 ${SQRT} ${PLACE}"
  fi
  echo "MMchecker ${Lsms} ${Rsms} ${Psms} ${MMFLAGS} "

  local mkn=`MMchecker ${Lsms} ${Rsms} ${Psms} ${MMFLAGS} |& grep '#'`
  echo $mkn
  if [[ "$mkn" == *"ERROR"* ]]; then
    exit 1;
  fi
else
  echo "<$m;$k;$n> algorithm of rank $r."
fi

Lmat=`dirname $Lsms`/`basename $Lsms .sms`
Rmat=`dirname $Lsms`/`basename $Rsms .sms`
Pmat=`dirname $Lsms`/`basename $Psms .sms`

#echo "$Lmat $Rmat $Pmat"
Lslp=${Lmat}.slp
Rslp=${Rmat}.slp
Pslp=${Pmat}.slp

##########
# Computation of the straight-line programs

for mat in $Lmat $Rmat
do
    Msms=${mat}.sms
    Mslp=${mat}.slp
    echo "Generating ${Mslp} with flags: ${OPTFLAGS}"
    optimizer ${OPTFLAGS} ${Msms} | compacter -s > ${Mslp}
done

echo "Generating ${Pslp}, by transposition, with flags: ${OPTFLAGS}"
matrix-transpose ${Psms} | optimizer ${OPTFLAGS} | transpozer | compacter -s > ${Pslp}


for mat in $Lmat $Rmat $Pmat
do
    Msms=${mat}.sms
    Mslp=${mat}.slp
    pmc=`PMchecker -M ${Msms} ${Mslp} |& grep '#'`
    echo $pmc
    if [[ "$pmc" == *"ERROR"* ]]; then
	exit 1;
    fi
done
}

##########
# Do generate the SLPs:

sms2slp ${Lsms} ${Rsms} ${Psms} ${MMCHECK} ${SQRT} ${PLACE}

##########
# Function replacing the placeholder by sqrt, with precomputed constants
function PlaceHolder {
    local SQRT=$1
    local PLACE=$2
    shift 2
    if [[ "$SQRT" -ne 0 ]]; then
	DBLP=$(( ${PLACE} * 2 ))
	sed -i "s/${DBLP}/2\*${PLACE}/g;s/${PLACE}/sqrt(${SQRT})/g;s/sqrt(${SQRT})\*2\/3/SQRT${SQRT}f2o3/g;s/2\*sqrt(${SQRT})\/3/SQRT${SQRT}f2o3/g;s/sqrt(${SQRT})\/2/SQRT${SQRT}o2/g;s/sqrt(${SQRT})\/3/SQRT${SQRT}o3/g;s/function .*/&\nSQRT${SQRT}o2=sqrt(${SQRT})\/2;\nSQRT${SQRT}o3=sqrt(${SQRT})\/3;\nSQRT${SQRT}f2o3=sqrt(${SQRT})\*2\/3;\n/" $*
    fi
}

##########
# Gathering of all straight-line program/variables to produce a Matlab program:

if [[ "$CoBTYPE" -eq 1 ]]; then
    echo "Generating ${File}}_CoBL.m ${File}_CoBR.m ${File}_ICoB.m change of bases:"
    ./CoB.rpl ${Lslp} ${m} ${k} ${File}_CoBL
    ./CoB.rpl ${Rslp} ${k} ${n} ${File}_CoBR
    ./CoB.rpl ${Pslp} ${m} ${n} ${File}_ICoB

    PlaceHolder ${SQRT} ${PLACE} ${File}_CoBL.m ${File}_CoBR.m ${File}_ICoB.m
else
    echo "Generating ${File}.m with ${m}x${k}x${n} of rank ${r}:"
    ./MM.rpl ${Lslp} ${Rslp} ${Pslp} ${m} ${k} ${n} ${r} ${File}

    PlaceHolder ${SQRT} ${PLACE} ${File}_${m}_${k}_${n}.m
fi
