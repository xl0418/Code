CFLAGS = -std=c++11 -ffast-math -march=native -O2 -Wall
OBJS = gdm.o sampler.o exp_expr.o iomatlab.o jc.o main.o

all: ${OBJS}
	g++ -o jc ${CFLAGS} ${OBJS}

gdm.o: gdm.cpp gdm.h
	g++ ${CFLAGS} -c gdm.cpp

sampler.o: sampler.cpp sampler.h
	g++ ${CFLAGS} -c sampler.cpp

exp_expr.o: exp_expr.cpp exp_expr.h
	g++ ${CFLAGS} -c exp_expr.cpp

iomatlab.o: iomatlab.cpp iomatlab.h
	g++ ${CFLAGS} -c iomatlab.cpp

jc.o: jc.cpp jc.h
	g++ ${CFLAGS} -c jc.cpp

main.o: main.cpp
	g++ ${CFLAGS} -c main.cpp

clean:
	rm -f jc ${OBJS}


