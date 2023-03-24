run: test
	./main
test:
	cc math.c -O3 -DMATHC_RUN -lm -o main
clean: 
	rm main
