CFLAGS=--ansi --pedantic -O2 -std=c++11 
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../src/
LDFLAGS=  ../src/*.o

all: epifire gsl ex8_event_driven_SEIRS

epifire: 
	$(MAKE) -C ../src/

gsl:
	$(MAKE) -C ../gsl_subset/

NUCOVID: ../src/Event_Driven_NUCOVID.h ../src/Utility.h NUCOVID.cpp
	g++ $(CFLAGS) NUCOVID.cpp $(INCLUDE) $(LDFLAGS) -o NUCOVID

model: ../src/Event_Driven_NUCOVID.h ../src/Utility.h ../src/Time_Series.h NUCOVID_2agegrp_modwave.cpp
	g++ $(CFLAGS) NUCOVID_2agegrp_modwave.cpp $(INCLUDE) $(LDFLAGS) -o model

model1g: ../src/Event_Driven_NUCOVID2.h ../src/Utility.h ../src/Time_Series.h NUCOVID_single_modwave.cpp
	g++ $(CFLAGS) NUCOVID_single_modwave.cpp $(INCLUDE) $(LDFLAGS) -o model1g

model2g: ../src/Event_Driven_NUCOVID2.h ../src/Utility.h ../src/Time_Series.h NUCOVID2_2agegrp_modwave.cpp
	g++ $(CFLAGS) NUCOVID2_2agegrp_modwave.cpp $(INCLUDE) $(LDFLAGS) -o model2g

clean:
	rm -f test_network ex1_mass_action ex2_percolation ex3_chain_binomial ex4_dynamic_net ex5_diff_eq ex6_network_diff_eq ex7_gillespie_network_SEIRS
