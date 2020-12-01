bicycle:func.cpp main.cpp
	#g++ -std=c++11 `itpp-config --cflags` -o bicycle.out func.cpp main.cpp  `itpp-config --libs`
	g++ -std=c++11 `itpp-config --cflags` -pthread -march=native -o bicycle.out func.cpp main.cpp  `itpp-config --libs` -lntl -lgmp -lm 
clean:
	rm *.out