CFLAGS = -std=c++11 -ffast-math -march=native -O2 -Wall
OBJS = gdm.o exp_expr.o iomatlab.o jc.o jc_neutral.o jc_implicit.o neutral.o main.o

all: ${OBJS}
	g++ -o jcimp ${CFLAGS} ${OBJS}

gdm.o:
	g++ ${CFLAGS} -c gdm.cpp

exp_expr.o: exp_expr.cpp exp_expr.h
	g++ ${CFLAGS} -c exp_expr.cpp

iomatlab.o: iomatlab.cpp iomatlab.h
	g++ ${CFLAGS} -c iomatlab.cpp

jc.o: jc.cpp jc.h
	g++ ${CFLAGS} -c jc.cpp

jc_implicit.o: jc_implicit.cpp
	g++ ${CFLAGS} -c jc_implicit.cpp

jc_neutral.o: jc_neutral.cpp
	g++ ${CFLAGS} -c jc_neutral.cpp

neutral.o: neutral.cpp neutral.h
	g++ ${CFLAGS} -c neutral.cpp

main.o: main.cpp
	g++ ${CFLAGS} -c main.cpp

clean:
	rm -f jcimp ${OBJS}


