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



##########
# L,R,P matrices
#
Lsms=$1
Rsms=$2
Psms=$3
File=$4

Lmat=`dirname $Lsms`/`basename $Lsms .sms`
Rmat=`dirname $Lsms`/`basename $Rsms .sms`
Pmat=`dirname $Lsms`/`basename $Psms .sms`
#echo "$Lmat $Rmat $Pmat"


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
# ==========================================================================


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
##########
# Function: Matrix Multiplication Check
#
function MMcheck {
    local Lsms=$1
    local Rsms=$2
    local Psms=$3
    local SQRT=$4
    local PLACE=$5

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
}
# ==========================================================================

# ==========================================================================
##########
# Function: Generator of SLPs from 3 (2 direct and 1 transposed) sms matrices
#
function sms2slp {
    local Lmat=$1
    local Rmat=$2
    local Pmat=$3
    local MMCHECK=$4
    local SQRT=$5
    local PLACE=$6

    local Lsms=${Lmat}.sms
    local Rsms=${Rmat}.sms
    local Psms=${Pmat}.sms

    local Lslp=${Lmat}.slp
    local Rslp=${Rmat}.slp
    local Pslp=${Pmat}.slp

    ##########
    # Do represent a matrix multiplication algorithm
    #
    if [[ "$MMCHECK" -eq 1 ]]; then
	MMcheck ${Lsms} ${Rsms} ${Psms} ${SQRT} ${PLACE}
    else
	echo "<$m;$k;$n> algorithm of rank $r."
    fi

    ##########
    # Computation of the straight-line programs
    #
    for mat in $Lmat $Rmat
      do
      local Msms=${mat}.sms
      local Mslp=${mat}.slp
      echo "Generating ${Mslp} with flags: ${OPTFLAGS}"
      optimizer ${OPTFLAGS} ${Msms} | compacter -s > ${Mslp}
    done

    echo "Generating ${Pslp}, by transposition, with flags: ${OPTFLAGS}"
    matrix-transpose ${Psms} | optimizer ${OPTFLAGS} | transpozer | compacter -s > ${Pslp}


    ##########
    # Verifications
    #
    for mat in $Lmat $Rmat $Pmat
      do
      local Msms=${mat}.sms
      local Mslp=${mat}.slp
      local pmc=`PMchecker -M ${Msms} ${Mslp} |& grep '#'`
      echo $pmc
      if [[ "$pmc" == *"ERROR"* ]]; then
	  exit 1;
      fi
    done
}
# ==========================================================================

# ==========================================================================
##########
# Function: replacing the placeholder by sqrt, with precomputed constants
#
function PlaceHolder {
    local SQRT=$1
    local PLACE=$2
    shift 2
    if [[ "$SQRT" -ne 0 ]]; then
	DBLP=$(( ${PLACE} * 2 ))
	sed -i "s/${DBLP}/2\*${PLACE}/g;s/${PLACE}/sqrt(${SQRT})/g;s/sqrt(${SQRT})\*2\/3/SQRT${SQRT}f2o3/g;s/2\*sqrt(${SQRT})\/3/SQRT${SQRT}f2o3/g;s/sqrt(${SQRT})\/2/SQRT${SQRT}o2/g;s/sqrt(${SQRT})\/3/SQRT${SQRT}o3/g;s/function .*/&\nSQRT${SQRT}o2=sqrt(${SQRT})\/2;\nSQRT${SQRT}o3=sqrt(${SQRT})\/3;\nSQRT${SQRT}f2o3=sqrt(${SQRT})\*2\/3;\n/" $*
    fi
}
# ==========================================================================




# ==========================================================================
##########
# Function: 
#
function slp2MMmpl {
    local Lslp=$1
    local Rslp=$2
    local Pslp=$3
    local m=$4
    local k=$5
    local n=$6
    local r=$7
    local suffix=$8

    filename="${m}x${k}x${n}_${r}_${suffix}.mpl"

    echo "# Left SLP" > ${filename}

echo "    (./replacer ${Lslp} -M i o A 3 3 40 1) >> ${filename} "
    (./replacer ${Lslp} -M i o A ${m} ${k} ${r} 1) >> ${filename}

    echo "# Right SLP" >> ${filename}
    (./replacer ${Rslp} -M i o B ${k} ${n} ${r} 1) >> ${filename}
    
    rmun=$((r-1))
    echo "# Inner products: ${rmun}" >> ${filename}
    for i in $(seq 0 ${rmun})
      do
      (echo -n "iC$i:=oA$i * oB$i; ") >> ${filename}
    done
    echo  >> ${filename}
    
    echo "# Post SLP" >> ${filename}
    (./replacer ${Pslp} -M i o C ${m} ${n} ${r} 0) >> ${filename}
    
    echo "# Check" >> ${filename}
    echo "Errors := LinearAlgebra:-Map(expand, C - ((Matrix(${m}, ${k}, symbol = 'A')) . (Matrix(${k}, ${n}, symbol = 'B')))); NumErrors:=LinearAlgebra:-Rank(Errors);" >> ${filename}

    if [ "${MAPLEHERE}" = true ] ; then
	tmpfile=/tmp/fdt_s2m.$$
	cat ${filename} | ${MAPLEPROG} | tee ${tmpfile}
	success=`grep "NumErrors := 0" ${tmpfile} | wc -l`
	if [[ "${success}" -ne 1 ]]; then
	    echo -e "\033[1;31m**** ERRORS in: cat ${filename} | ${MAPLEPROG} ***\033[0m"
	else
	   echo -e "\033[1;32mSUCCESS: ${filename} is a MM algorithm.\033[0m"
	fi
    fi
}




MMcheck ${Lsms} ${Rsms} ${Psms} ${SQRT} ${PLACE}

BUILDSLP=0

if [[ "$BUILDSLP" -eq 1 ]]; then
	# Do generate the SLPs:
    sms2slp ${Lmat} ${Rmat} ${Pmat} ${MMCHECK} ${SQRT} ${PLACE}
fi

Lslp=${Lmat}.slp
Rslp=${Rmat}.slp
Pslp=${Pmat}.slp


	# Bilinear algorithm for matrix multiplication only
slp2MMmpl ${Lslp} ${Rslp} ${Pslp} ${m} ${k} ${n} ${r} ${File}
