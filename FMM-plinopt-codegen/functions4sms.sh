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
# Function: FMM-codegen check
#
function FMMcodegenPresent {
    for prg in replacer CoB.rpl MM.rpl
      do
      if ! command -v `dirname $0`/${prg} &> /dev/null
      then
	  echo -e "\033[1;31m**** ERROR\033[0m FMM code generators (${prg}) could not be found. \033[1;31m***\033[0m"
	  exit 1;
      fi
    done
    echo "# FMM-codegen is up and running"
}


# ==========================================================================
##########
# Function: PLinOpt check
#
function PLinOptPresent {
    for prg in MMchecker optimizer matrix-transpose SLPchecker factorizer sms2pretty columns-swap
      do
      if ! command -v ${prg} &> /dev/null
      then
	  echo -e "\033[1;31m**** ERROR\033[0m plinopt/bin executables (${prg}) could not be found. \033[1;31m***\033[0m"
	  exit 1;
      fi
    done
    echo "# PLinOpt is up and running"
}

# ==========================================================================
##########
# Function: Matrix Multiplication Check
#
function MMcheck {
    local Lsms=$1
    local Rsms=$2
    local Psms=$3
    local REPL=$4
    local EXPO=$5
    local SQRT=$6

    local MMFLAGS=""
    if [[ "${SQRT}" -ne 0 ]]; then
	local MMFLAGS="-b 128 -r ${REPL} ${EXPO} ${SQRT}"
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
    local REPL=$5
    local EXPO=$6
    local SQRT=$7
    local OPTFLAGS=$8

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
	MMcheck ${Lsms} ${Rsms} ${Psms} ${REPL} ${EXPO} ${SQRT}
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
    local REPL=$1
    local EXPO=$2
    local SQRT=$3
    shift 3
    if [[ "${SQRT}" -ne 0 ]]; then
	DBLP=$(( ${REPL} * 2 ))
	sed -i "s/${DBLP}/2\*${REPL}/g;s/${REPL}/nthroot(${SQRT},${EXPO})/g;s/nthroot(${SQRT},${EXPO})\*2\/3/NTH${EXPO}ROOT${SQRT}f2o3/g;s/2\*nthroot(${SQRT},${EXPO})\/3/NTH${EXPO}ROOT${SQRT}f2o3/g;s/1\/nthroot(${SQRT},${EXPO})/NTH${EXPO}ROOT${SQRT}t1o/g;s/2\/nthroot(${SQRT},${EXPO})/NTH${EXPO}ROOT${SQRT}t2o/g;s/nthroot(${SQRT},${EXPO})\/2/NTH${EXPO}ROOT${SQRT}o2/g;s/nthroot(${SQRT},${EXPO})\/4/NTH${EXPO}ROOT${SQRT}o4/g;s/nthroot(${SQRT},${EXPO})\/8/NTH${EXPO}ROOT${SQRT}o8/g;s/nthroot(${SQRT},${EXPO})\/16/NTH${EXPO}ROOT${SQRT}o16/g;s/nthroot(${SQRT},${EXPO})\/3/NTH${EXPO}ROOT${SQRT}o3/g;s/nthroot(${SQRT},${EXPO})/NTH${EXPO}ROOT${SQRT}once/g;s/function .*/&\nNTH${EXPO}ROOT${SQRT}once=nthroot(${SQRT},${EXPO});\nNTH${EXPO}ROOT${SQRT}o2=nthroot(${SQRT},${EXPO})\/2;\nNTH${EXPO}ROOT${SQRT}t2o=2\/nthroot(${SQRT},${EXPO});\nNTH${EXPO}ROOT${SQRT}t1o=1\/nthroot(${SQRT},${EXPO});\nNTH${EXPO}ROOT${SQRT}o4=nthroot(${SQRT},${EXPO})\/4;\nNTH${EXPO}ROOT${SQRT}o8=nthroot(${SQRT},${EXPO})\/8;\nNTH${EXPO}ROOT${SQRT}o16=nthroot(${SQRT},${EXPO})\/16;\nNTH${EXPO}ROOT${SQRT}o3=nthroot(${SQRT},${EXPO})\/3;\nNTH${EXPO}ROOT${SQRT}f2o3=nthroot(${SQRT},${EXPO})\*2\/3;\n/" $*
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
    local REPL=$8
    local EXPO=$9
    local SQRT=${10}
    local File=${11}
    echo "# Generating ${File}_CoBL.m ${File}_CoBR.m ${File}_ICoB.m change of bases:"
    `dirname $0`/CoB.rpl ${Lslp} ${m} ${k} ${File}_CoBL
    `dirname $0`/CoB.rpl ${Rslp} ${k} ${n} ${File}_CoBR
    `dirname $0`/CoB.rpl ${Pslp} ${m} ${n} ${File}_ICoB

    PlaceHolder ${REPL} ${EXPO} ${SQRT} ${File}_CoBL.m ${File}_CoBR.m ${File}_ICoB.m
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
    local REPL=$8
    local EXPO=$9
    local SQRT=${10}
    local File=${11}
    echo "# Generating ${File}.m with ${m}x${k}x${n} of rank ${r}:"
    `dirname $0`/MM.rpl ${Lslp} ${Rslp} ${Pslp} ${m} ${k} ${n} ${r} ${File} 1

echo "#### PlaceHolder ${REPL} ${EXPO} ${SQRT} ${File}_${m}_${k}_${n}.m ####"
    PlaceHolder ${REPL} ${EXPO} ${SQRT} ${File}_${m}_${k}_${n}.m
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
    local REPL=$9
    local EXPO=${10}
    local SQRT=${11}
    MODULO=$(( ${REPL} ** ${EXPO} - ${SQRT} ))

    filename="${m}x${k}x${n}_${r}_${suffix}.mpl"

    echo "# Left SLP" > ${filename}
    (`dirname $0`/replacer ${Lslp} -M i o A ${m} ${k} ${r} 1 | sed 's/b/w/g;s/x/s/g;s/t/x/g;s/v/t/g;s/z/v/g;s/g/u/g;s/oA/l/g') >> ${filename}

    echo "# Right SLP" >> ${filename}
    (`dirname $0`/replacer ${Rslp} -M i o B ${k} ${n} ${r} 1| sed 's/t/y/g;s/b/g/g;s/x/c/g;s/v/d/g;s/g/e/g;s/z/f/g;s/r/b/g;s/oB/r/g') >> ${filename}

    rmun=$((r-1))
    echo "# Inner products: ${rmun}" >> ${filename}
    for i in $(seq 0 ${rmun})
      do
      (echo -n "p${i}:=l${i}*r${i}; ") >> ${filename}
    done
    echo -e '\n'  >> ${filename}

    echo "C:=Matrix(${m},${n}):" >> ${filename}
    echo -e '\n# Post SLP' >> ${filename}
    (`dirname $0`/replacer ${Pslp} -M i o C ${m} ${n} ${r} 0| sed 's/z/n/g;s/r/h/g;s/x/j/g;s/v/k/g;s/g/m/g;s/t/z/g;s/b/q/g;s/iC/p/g') >> ${filename}

    echo "# Check" >> ${filename}
    modcomp="";
    if [[ "${SQRT}" -ne 0 ]]; then
	modcomp=" mod ${MODULO}";
    fi

    echo "Errors := LinearAlgebra:-Map(expand, C - ((Matrix(${m}, ${k}, symbol = 'A')) . (Matrix(${k}, ${n}, symbol = 'B'))))"${modcomp}"; NumErrors:=LinearAlgebra:-Rank(Errors);" >> ${filename}

    if [ "${MAPLEHERE}" = true ] ; then
	tmpfile=/tmp/fdt_s2m.$$
	cat ${filename} | ${MAPLEPROG} | tee ${tmpfile} | tail
	success=`grep "NumErrors := 0" ${tmpfile} | wc -l`
	if [[ "${success}" -ne 1 ]]; then
	    cat ${tmpfile}
	    echo -e "\033[1;31m**** ERRORS in: cat ${filename} | ${MAPLEPROG} ***\033[0m"
	    exit 1;
	else
	    echo -en "\033[1;32mSUCCESS:\033[0m "
	fi
	\rm ${tmpfile}
    fi

    echo "${filename} is a MM algorithm."
}
# ==========================================================================
