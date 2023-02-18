run: main_run
	./main
main_lib:
	cc main.c -O3 -c -o main
main_run:
	cc main.c -O3 -DMATH_RUN -lm -o main
clean: 
	rm main
