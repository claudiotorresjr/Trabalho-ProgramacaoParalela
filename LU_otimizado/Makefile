		PROG   = sisLinear
		OBJS   = $(PROG).o fatLU.o

#		LIKWID = /home/soft/likwid
#	LIKWID_FLAGS = -DLIKWID_PERFMON -I$(LIKWID)/include
#	LIKWID_LIBS = -L$(LIKWID)/lib -llikwid

		CC = gcc -std=c11 -O3 -mavx -march=native -fopenmp
	CFLAGS = -Wall #$(LIKWID_FLAGS)
	LFLAGS = -lm   #$(LIKWID_LIBS) 

.PHONY: clean limpa purge faxina distclean debug avx

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak *.tmp

purge distclean:   clean
	@rm -f  $(PROG) *.o core a.out
