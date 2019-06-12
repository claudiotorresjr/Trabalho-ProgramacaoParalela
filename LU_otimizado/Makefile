		PROG   = sisLinear
		OBJS   = $(PROG).o bloco_certo.o  #fatLU.o #avx_bloco.o

#		LIKWID = /home/soft/likwid
#	LIKWID_FLAGS = -DLIKWID_PERFMON -I$(LIKWID)/include
#	LIKWID_LIBS = -L$(LIKWID)/lib -llikwid

		CC = gcc -O3 -mavx -fopenmp
	CFLAGS = -Wall #$(LIKWID_FLAGS)
	LFLAGS = #-lm   $(LIKWID_LIBS) 

.PHONY: clean limpa purge faxina distclean debug avx

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak *.tmp

purge distclean:   clean
	@rm -f  $(PROG) *.o core a.out
