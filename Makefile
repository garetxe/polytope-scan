PALPDIR = /home/inaki/cvs/palp
BASEDIR = /home/inaki/polytopes/4d
LIBSDIR = /home/inaki/cvs/sage/local/lib
INCLUDEDIR = /home/inaki/cvs/sage/local/include/python2.6

scan: scan.c
	gcc -oscan scan.c -Wall -Wextra -O2 -g -lpthread -L$(LIBSDIR) -I$(INCLUDEDIR) $(LIBSDIR)/libcddgmp.a -lgmp

test: scan
	time ./scan 2 "$(PALPDIR)/class.x -b2a -di $(BASEDIR)/zzdb" "$(PALPDIR)/poly.x"

all: scan
