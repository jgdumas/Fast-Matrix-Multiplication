CXXFLAGS += ${OPTFLAGS} -O3
SHELL=/bin/bash

all: replacer

check: mtlcheck mplcheck

mtlcheck: sms2matlab.sh
	echo n | ./sms2matlab.sh data/{L,R,P}o.sms data/DPS -p 3 1013

mplcheck: sms2maple.sh
	echo n | ./sms2maple.sh data/{L,R,P}_3x3x6.sms check
	echo n | ./sms2maple.sh data/4x4x4_48_rational_{L,R,P}.sms check
	echo n | ./sms2maple.sh data/3x4x7_63_rational_{L,R,P}.sms check
