CC=mpicc -O3 -march=native -g
OBJECTS=main.o $(patsubst src/%.c, obj/%.o, $(wildcard src/*.c))

main.x: $(OBJECTS)
	$(CC) $^ -o $@ -lm

main.o: main.c
	$(CC) -c $^

obj/%.o: src/%.c
	$(CC) -c $^ -o $@
clean:
	rm -rf *.o *.x
	rm -rf obj/*.o
