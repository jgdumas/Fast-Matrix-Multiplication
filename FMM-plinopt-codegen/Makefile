CXXFLAGS += ${OPTFLAGS} -O3

all: replacer

check: sms2maple.sh
	./sms2maple.sh data/L_3x3x6.sms data/R_3x3x6.sms data/P_3x3x6.sms check
	./sms2maple.sh data/4x4x4_48_rational_L.sms data/4x4x4_48_rational_R.sms data/4x4x4_48_rational_P.sms check
	./sms2maple.sh data/3x4x7_63_rational_L.sms data/3x4x7_63_rational_R.sms data/3x4x7_63_rational_P.sms check