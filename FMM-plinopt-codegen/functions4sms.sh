#!/bin/bash
# ==========================================================================
# Fast-Matrix-Multiplication library
#   https://github.com/jgdumas/Fast-Matrix-Multiplication
#   Accurate fast matrix multiplications via recursion
# ==========================================================================
# function4sms: tools for matlab/maple program generation
# ==========================================================================
#   Authors:
#   [J-G. Dumas, C. Pernet, A. Sedoglavic;
#    Strassen's algorithm is not optimally accurate;
#    ISSAC 2024, Raleigh, NC USA, pp. 254-263.
#    https://hal.science/hal-04441653]
# ==========================================================================

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
    echo "# MMchecker ${Lsms} ${Rsms} ${Psms} ${MMFLAGS}"
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
    local OPTFLAGS=$7

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
	echo "# <$m;$k;$n> algorithm of rank $r."
    fi

    ##########
    # Computation of the straight-line programs
    #
    BSLP=0

    if [ -f ${Lslp} ] || [ -f ${Rslp} ] || [ -f ${Lslp} ]; then
	read -p "# Overwrite ${Lslp}, ${Rslp} and ${Pslp}? [y/N] " RESP
	if [[ "${RESP}" == "y" ]]; then
	    BSLP=1
	fi
    else
	BSLP=1
    fi

    if [[ "$BSLP" -eq 1 ]]; then
	# Do generate the SLPs:
	for mat in $Lmat $Rmat
	  do
	  local Msms=${mat}.sms
	  local Mslp=${mat}.slp
	  echo "# Generating ${Mslp} with flags: ${OPTFLAGS}"
	  optimizer ${OPTFLAGS} ${Msms} | compacter -s > ${Mslp}
	done

	echo "# Generating ${Pslp}, by transposition, with flags: ${OPTFLAGS}"
	matrix-transpose ${Psms} | optimizer ${OPTFLAGS} | transpozer | compacter -s > ${Pslp}
    fi

    ##########
    # Verifications
    #
    for mat in $Lmat $Rmat $Pmat
      do
      local Msms=${mat}.sms
      local Mslp=${mat}.slp
      local pmc=`SLPchecker -M ${Msms} ${Mslp} |& grep '#'`
      echo $pmc | sed 's/#[^#]/\n&/g'
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
# Function: Gathering of all CoB SLPs to produce a Matlab program
#
function slp2CBm {
    local Lslp=$1
    local Rslp=$2
    local Pslp=$3
    local m=$4
    local k=$5
    local n=$6
    local r=$7
    local File=$8
    echo "# Generating ${File}_CoBL.m ${File}_CoBR.m ${File}_ICoB.m change of bases:"
    ./CoB.rpl ${Lslp} ${m} ${k} ${File}_CoBL
    ./CoB.rpl ${Rslp} ${k} ${n} ${File}_CoBR
    ./CoB.rpl ${Pslp} ${m} ${n} ${File}_ICoB

    PlaceHolder ${SQRT} ${PLACE} ${File}_CoBL.m ${File}_CoBR.m ${File}_ICoB.m
}
# ==========================================================================


# ==========================================================================
##########
# Function: Gathering of SLPs to form a Matlab program
#
function slp2MMm {
    local Lslp=$1
    local Rslp=$2
    local Pslp=$3
    local m=$4
    local k=$5
    local n=$6
    local r=$7
    local File=$8
    echo "# Generating ${File}.m with ${m}x${k}x${n} of rank ${r}:"
    ./MM.rpl ${Lslp} ${Rslp} ${Pslp} ${m} ${k} ${n} ${r} ${File} 1
    PlaceHolder ${SQRT} ${PLACE} ${File}_${m}_${k}_${n}.m
}
# ==========================================================================


# ==========================================================================
##########
# Function: combine two squential programs to check with linear operator
#
function combPMcheck {
    local mat=$1
    local pri=$2
    local pro=$3
    local usedvar=`sed 's/[^[:alpha:]]//g' ${pri} ${pro} | grep -o . | sort -u | tr -d "\n"`
    local freechar=`echo "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" | sed "s/${usedvar}//g" | grep -o . | head -1`
    sed -r "s/(o)([0-9]*)/${freechar}\2/g" ${pri} > s2m.slp
    sed -r "s/(i)([0-9]*)/${freechar}\2/g" ${pro} >> s2m.slp

    local pmc=`SLPchecker -M ${mat} s2m.slp |& grep '#'`
    echo $pmc
    if [[ "$pmc" == *"ERROR"* ]]; then
	echo "# Combined ${pri} ${pro}: SLPchecker -M ${mat} s2m.slp"
	exit 1;
    else
	\rm s2m.slp
    fi
}
# ==========================================================================



# ==========================================================================
##########
# Function: Gathering of SLPs to form a Maple program
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
    (./replacer ${Lslp} -M i o A ${m} ${k} ${r} 1 | sed 's/b/w/g;s/x/s/g;s/t/x/g;s/v/t/g;s/z/v/g;s/g/u/g;s/oA/l/g') >> ${filename}

    echo "# Right SLP" >> ${filename}
    (./replacer ${Rslp} -M i o B ${k} ${n} ${r} 1| sed 's/t/y/g;s/b/g/g;s/x/c/g;s/v/d/g;s/g/e/g;s/z/f/g;s/r/b/g;s/oB/r/g') >> ${filename}

    rmun=$((r-1))
    echo "# Inner products: ${rmun}" >> ${filename}
    for i in $(seq 0 ${rmun})
      do
      (echo -n "p${i}:=l${i}*r${i}; ") >> ${filename}
    done
    echo -e '\n'  >> ${filename}

    echo "C:=Matrix(${m},${n}):" >> ${filename}
    echo -e '\n# Post SLP' >> ${filename}
    (./replacer ${Pslp} -M i o C ${m} ${n} ${r} 0| sed 's/z/n/g;s/r/h/g;s/x/j/g;s/v/k/g;s/g/m/g;s/t/z/g;s/b/q/g;s/iC/p/g') >> ${filename}

    echo "# Check" >> ${filename}
    echo "Errors := LinearAlgebra:-Map(expand, C - ((Matrix(${m}, ${k}, symbol = 'A')) . (Matrix(${k}, ${n}, symbol = 'B')))); NumErrors:=LinearAlgebra:-Rank(Errors);" >> ${filename}

    if [ "${MAPLEHERE}" = true ] ; then
	tmpfile=/tmp/fdt_s2m.$$
	cat ${filename} | ${MAPLEPROG} | tee ${tmpfile} | tail
	success=`grep "NumErrors := 0" ${tmpfile} | wc -l`
	if [[ "${success}" -ne 1 ]]; then
	    echo -e "\033[1;31m**** ERRORS in: cat ${filename} | ${MAPLEPROG} ***\033[0m"
	    exit 1;
	else
	    echo -en "\033[1;32mSUCCESS:\033[0m "
	fi
    fi

    echo "${filename} is a MM algorithm."
}
# ==========================================================================
