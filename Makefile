run: main_run
	./main
main:
	cc main.c -O3 -lm -o main
main_run:
	cc main.c -O3 -DMATH_RUN -lm -o main
clean: 
	rm main
