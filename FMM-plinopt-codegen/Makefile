CXXFLAGS += ${OPTFLAGS} -O3
SHELL=/bin/bash
BIN=replacer
SCR=MM.rpl CoB.rpl

all: ${BIN}

check: mtlcheck mplcheck

mtlcheck: sms2matlab.sh ${SCR} ${BIN}
	echo n | ./sms2matlab.sh data/2x2x2_7_DPS-accurate_{L,R,P}.sms data/DPS -p 3 1013

mplcheck: sms2maple.sh ${SCR} ${BIN}
	echo n | ./sms2maple.sh data/3x3x6_40_DPS-accurate_{L,R,P}.sms check
	echo n | ./sms2maple.sh data/4x4x4_48_rational_{L,R,P}.sms check
	echo n | ./sms2maple.sh data/3x4x7_63_rational_{L,R,P}.sms check

clean:
	- \rm ${BIN}

