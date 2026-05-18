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
OPTFLAGS="-O 10"
FCTFLAGS=""
MMCHECK=1
CoBTYPE=0
REPL=0
EXPO=0
SQRT=0
MODU=0
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
      echo "Usage: $0 [-O #|-r # # #|-m|-n|-c|-a|-f #] L.sms R.sms P.sms name"
      echo "  generates matlab program name.m from L,R,P matrices."
      echo "  -m: L,R,P are a matrix multiplication algorithm (default)."
      echo "  -c: L,R,P are change of bases matrices (default MM)."
      echo "  -a: L,R,P are first factorized with alternative bases."
      echo "  -f k: column dimension of factorized bases (default is column dim.)."
      echo "  -m/-n: L,R,P are checked/not checked as a mat. mul. (default no)."
      echo "  -O N: optimizer with N loops (default is ${OPTFLAGS})."
      echo "  -q q: compute modulo q (default none)."
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


LAmat="`dirname $Lsms`/`basename $Lsms _L.sms`-ALT_L"
RAmat="`dirname $Rsms`/`basename $Rsms _R.sms`-ALT_R"
PAmat="`dirname $Psms`/`basename $Psms _P.sms`-ALT_P"
LCmat="`dirname $Lsms`/`basename $Lsms _L.sms`-CoB_L"
RCmat="`dirname $Rsms`/`basename $Rsms _R.sms`-CoB_R"
PCmat="`dirname $Psms`/`basename $Psms _P.sms`-CoB_P"
LAsms="${LAmat}.sms"
RAsms="${RAmat}.sms"
PAsms="${PAmat}.sms"
LCsms="${LCmat}.sms"
RCsms="${RCmat}.sms"
PCsms="${PCmat}.sms"
LAslp="${LAmat}.slp"
RAslp="${RAmat}.slp"
PAslp="${PAmat}.slp"
LCslp="${LCmat}.slp"
RCslp="${RCmat}.slp"
PCslp="${PCmat}.slp"

# echo "### $LAsms $RAsms $PAsms"
# echo "### $LCsms $RCsms $PCsms"
# echo "### $LAslp $RAslp $PAslp"
# echo "### $LCslp $RCslp $PCslp"

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

    MMcheck ${Lsms} ${Rsms} ${Psms} ${REPL} ${EXPO} ${SQRT} ${MODU}

    OvwrCoB=0
    if [ -f ${LAsms} ] || [ -f ${RAsms} ] || [ -f ${PAsms} ] || [ -f ${LCsms} ] || [ -f ${RCsms} ] || [ -f ${PCsms} ]; then
	read -p "# Overwrite ${LAsms}, ${RAsms}, ${PAsms}, ${LCsms}, ${RCsms}, ${PCsms}? [y/N] " RESP
	if [[ "${RESP}" == "y" ]]; then
	    OvwrCoB=1
	fi
    else
	OvwrCoB=1
    fi

    if [[ "$OvwrCoB" -eq 1 ]]; then
	# Factor into: sparse x CoB
	echo "# Sparsifying into ${LAsms} ${RAsms} ${PAsms}; ${FCTFLAGS}"
	echo "#       with bases ${LCsms} ${RCsms} ${PCsms}."
	(factorizer -S ${OPTFLAGS} ${FCTFLAGS} ${Lsms} > ${LCsms}) |& grep -v '#' > ${LAsms}
	(factorizer -S ${OPTFLAGS} ${FCTFLAGS} ${Rsms} > ${RCsms}) |& grep -v '#' > ${RAsms}

	(matrix-transpose ${Psms} | factorizer -S ${OPTFLAGS} ${FCTFLAGS} > ${PCmat}_tC.sms) |& grep -v '#' > ${PCmat}_tA.sms
	matrix-transpose ${PCmat}_tC.sms > ${PCsms}
	matrix-transpose ${PCmat}_tA.sms > ${PAsms}
    else
	echo "# Using ${LAsms}, ${RAsms}, ${PAsms}, ${LCsms}, ${RCsms}, ${PCsms}"s
    fi

	# Do generate the SLPs:
    sms2slp ${LCmat} ${RCmat} ${PCmat} 0 ${REPL} ${EXPO} ${SQRT} ${MODU} "${OPTFLAGS}"
    sms2slp ${LAmat} ${RAmat} ${PAmat} 0 ${REPL} ${EXPO} ${SQRT} ${MODU} "${OPTFLAGS}"

	# Verifications
    combPMcheck ${Lsms} ${LCslp} ${LAslp}
    combPMcheck ${Rsms} ${RCslp} ${RAslp}
    combPMcheck ${Psms} ${PAslp} ${PCslp}

	# Produce the matlab programs from the SLPs, with CoB and bilinear algorithm

    slp2CBm ${LCslp} ${RCslp} ${PCslp} ${m} ${k} ${n} ${fl} ${fr} ${fp} ${REPL} ${EXPO} ${SQRT} ${File}
    if [[ "${fl}" -eq 0 ]]; then
	fl=$((m*k))
    fi
    if [[ "${fr}" -eq 0 ]]; then
	fr=$((k*n))
    fi
    if [[ "${fp}" -eq 0 ]]; then
	fp=$((m*n))
    fi
    slp2MMm ${LAslp} ${RAslp} ${PAslp} ${m} ${k} ${n} ${r} ${fl} ${fr} ${fp} ${REPL} ${EXPO} ${SQRT} ${File}_mul

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

echo -e "${GRE}SUCCESS: ${File}_alternative.m generated.${NC} "


else
    # Generic bilinear algorithm

	# Do generate the SLPs:
    sms2slp ${Lmat} ${Rmat} ${Pmat} ${MMCHECK} ${REPL} ${EXPO} ${SQRT} ${MODU} "${OPTFLAGS}"

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
