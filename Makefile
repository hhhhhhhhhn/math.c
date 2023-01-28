run: main
	./main
main:
	cc main.c -O3 -lm -o main
clean: 
	rm main
