CXXFLAGS += ${OPTFLAGS} -O3

all: replacer

check: mtlcheck mplcheck

mtlcheck: sms2matlab.sh
	./sms2matlab.sh data/Lo.sms data/Ro.sms data/Po.sms data/DPS -p 3 1013

mplcheck: sms2maple.sh
	echo n | ./sms2maple.sh data/L_3x3x6.sms data/R_3x3x6.sms data/P_3x3x6.sms check
	echo n | ./sms2maple.sh data/4x4x4_48_rational_L.sms data/4x4x4_48_rational_R.sms data/4x4x4_48_rational_P.sms check
	echo n | ./sms2maple.sh data/3x4x7_63_rational_L.sms data/3x4x7_63_rational_R.sms data/3x4x7_63_rational_P.sms check
