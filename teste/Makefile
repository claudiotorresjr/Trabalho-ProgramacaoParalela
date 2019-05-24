		PROG   = matriz_desempenho
		OBJS   = $(PROG).o matriz.o

		LIKWID = /home/soft/likwid
	LIKWID_FLAGS = -DLIKWID_PERFMON -I$(LIKWID)/include
	LIKWID_LIBS = -L$(LIKWID)/lib -llikwid

		CC = gcc -std=c11
	CFLAGS = $(LIKWID_FLAGS)
	LFLAGS = $(LIKWID_LIBS) -O3 -Wall -lm

.PHONY: clean limpa purge faxina distclean debug avx

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

avx: CFLAGS += -O3

clean:
	@rm -f *~ *.bak *.tmp

purge distclean:   clean
	@rm -f  $(PROG) *.o core a.out
