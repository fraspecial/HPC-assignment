CC=mpicc -O3 -march=native
CFLAGS= -c -I include
OBJECTS=main.o $(patsubst src/%.c, obj/%.o, $(wildcard src/*.c))

main.x: $(OBJECTS)
	$(CC) $^ -o $@

main.o: main.c
	$(CC) $(CFLAGS) $^
	obj/%.o: src/%.c
	$(CC) $(CFLAGS) $^ -o $@
clean:
	rm -rf *.o *.x
	rm -rf obj/*.o
