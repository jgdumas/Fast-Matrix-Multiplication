#!/bin/bash
# ==========================================================================
# Fast-Matrix-Multiplication library
#   https://github.com/jgdumas/Fast-Matrix-Multiplication
#   Accurate fast matrix multiplications via recursion
# ==========================================================================
# sms2matlab: generate matlab programs from an L,R,P bilinear algorithm
# ==========================================================================
#   Authors:
#   [J-G. Dumas, C. Pernet, A. Sedoglavic;
#    Strassen's algorithm is not optimally accurate;
#    ISSAC 2024, Raleigh, NC USA, pp. 254-263.
#    https://hal.science/hal-04441653]
# ==========================================================================


source `dirname $0`/functions4sms.sh

# ==========================================================================
##########
# Parsing args
#
MATS=()
OPTFLAGS="-E -N"
FCTFLAGS=""
MMCHECK=1
CoBTYPE=0
REPL=0
EXPO=0
SQRT=0
fl=0
fr=0
fp=0

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
    -c|-cob|--CoB|--ChangeofBasis)
      MMCHECK=0
      CoBTYPE=1
      shift # past argument
      ;;
    -n|-nm|-nmm|-nmmc|-nc|--NoMMcheck)
      MMCHECK=0
      shift # past argument
      ;;
    -a|-alt|--AlternativeBasis)
      ALTBASIS=1
      shift # past argument
      ;;
    -f|-k|-fact|--FactorizerFlags)
      fl=$2
      fr=$2
      fp=$2
      FCTFLAGS="-k $2"
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
      echo "Usage: $0 [-O #|-r # # #|-m|-n|-c|-a|-f #] L.sms R.sms P.sms name"
      echo "  generates matlab program name.m from L,R,P matrices."
      echo "  -m: L,R,P are a matrix multiplication algorithm (default)."
      echo "  -c: L,R,P are change of bases matrices (default MM)."
      echo "  -a: L,R,P are first factorized with alternative bases."
      echo "  -f k: column dimension of factorized bases (default is column dim.)."
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


# ==========================================================================
# ==========================================================================
##########
# Starting sms2matlab for -m/-n/-c options
#
if [[ "$ALTBASIS" -eq 1 ]]; then
    # Compute with alternative bases

    MMcheck ${Lsms} ${Rsms} ${Psms} ${REPL} ${EXPO} ${SQRT}

    OvwrCoB=0
    if [ -f ${Lmat}_A.sms ] || [ -f ${Rmat}_A.sms ] || [ -f ${Pmat}_A.sms ] || [ -f ${Lmat}_C.sms ] || [ -f ${Rmat}_C.sms ] || [ -f ${Pmat}_C.sms ]; then
	read -p "# Overwrite ${Lmat}_A.sms, ${Rmat}_A.sms, ${Pmat}_A.sms, ${Lmat}_C.sms, ${Rmat}_C.sms, ${Pmat}_C.sms? [y/N] " RESP
	if [[ "${RESP}" == "y" ]]; then
	    OvwrCoB=1
	fi
    else
	OvwrCoB=1
    fi

    if [[ "$OvwrCoB" -eq 1 ]]; then
	# Factor into: sparse x CoB
	echo "# Sparsifying into ${Lmat}_A.sms ${Rmat}_A.sms ${Pmat}_A.sms; ${FCTFLAGS}"
	echo "#       with bases ${Lmat}_C.sms ${Rmat}_C.sms ${Pmat}_C.sms."
	(factorizer -S ${FCTFLAGS} ${Lsms} > ${Lmat}_C.sms) |& grep -v '#' > ${Lmat}_A.sms
	(factorizer -S ${FCTFLAGS} ${Rsms} > ${Rmat}_C.sms) |& grep -v '#' > ${Rmat}_A.sms

	(matrix-transpose ${Psms} | factorizer -S ${FCTFLAGS} > ${Pmat}_tC.sms) |& grep -v '#' > ${Pmat}_tA.sms
	matrix-transpose ${Pmat}_tC.sms > ${Pmat}_C.sms
	matrix-transpose ${Pmat}_tA.sms > ${Pmat}_A.sms
    fi

	# Do generate the SLPs:
    sms2slp ${Lmat}_C ${Rmat}_C ${Pmat}_C 0 ${REPL} ${EXPO} ${SQRT} "${OPTFLAGS}"
    sms2slp ${Lmat}_A ${Rmat}_A ${Pmat}_A 0 ${REPL} ${EXPO} ${SQRT} "${OPTFLAGS}"

	# Verifications
    combPMcheck ${Lsms} ${Lmat}_C.slp ${Lmat}_A.slp
    combPMcheck ${Rsms} ${Rmat}_C.slp ${Rmat}_A.slp
    combPMcheck ${Psms} ${Pmat}_A.slp ${Pmat}_C.slp

	# Produce the matlab programs from the SLPs, with CoB and bilinear algorithm

    slp2CBm ${Lmat}_C.slp ${Rmat}_C.slp ${Pmat}_C.slp ${m} ${k} ${n} ${fl} ${fr} ${fp} ${REPL} ${EXPO} ${SQRT} ${File}
    if [[ "${fl}" -eq 0 ]]; then
	fl=$((m*k))
    fi
    if [[ "${fr}" -eq 0 ]]; then
	fr=$((k*n))
    fi
    if [[ "${fp}" -eq 0 ]]; then
	fp=$((m*n))
    fi
    slp2MMm ${Lmat}_A.slp ${Rmat}_A.slp ${Pmat}_A.slp ${m} ${k} ${n} ${r} ${fl} ${fr} ${fp} ${REPL} ${EXPO} ${SQRT} ${File}_mul

    Mfile=`basename ${File}`
    echo "# Generating alternative basis matlab program ${File}_alternative.m."

## CAT .......................
cat > ${File}_alternative.m<< EOF
function C = ${Mfile}_alternative(A, B, nmin, peeling, level)
%          Computes the product C = A*B, sparsified via alternative basis.
  if nargin < 3, nmin = 4; end    % Threshold to conventional
  if nargin < 4, peeling = 1; end % Static (1) or Dynamic (2) peeling
  if nargin < 5, level = 8; end   % Verbose level

  U = ${Mfile}_CoBL(A, nmin, peeling, level);   % Left change of basis
  V = ${Mfile}_CoBR(B, nmin, peeling, level);   % Right change of basis
  W = ${Mfile}_mul(U, V, nmin, peeling, level); % Sparse multiplication
  C = ${Mfile}_ICoB(W, nmin, peeling, level);   % Inverse change of basis
end
EOF
## ....................... EOF

echo -e "\033[1;32mSUCCESS: ${File}_alternative.m generated.\033[0m "


else
    # Generic bilinear algorithm

	# Do generate the SLPs:
    sms2slp ${Lmat} ${Rmat} ${Pmat} ${MMCHECK} ${REPL} ${EXPO} ${SQRT} "${OPTFLAGS}"

    Lslp=${Lmat}.slp
    Rslp=${Rmat}.slp
    Pslp=${Pmat}.slp

	# Produce the matlab program from the SLPs
    if [[ "$CoBTYPE" -eq 1 ]]; then
	# Change of bases only
	slp2CBm ${Lslp} ${Rslp} ${Pslp} ${m} ${k} ${n} ${fl} ${fr} ${fp} ${REPL} ${EXPO} ${SQRT} ${File}
    else
	# Bilinear algorithm for matrix multiplication only
	slp2MMm ${Lslp} ${Rslp} ${Pslp} ${m} ${k} ${n} ${r} ${fl} ${fr} ${fp} ${REPL} ${EXPO} ${SQRT} ${File}
    fi
fi

# ==========================================================================
# ==========================================================================
