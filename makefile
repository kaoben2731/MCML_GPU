MCML_GPU: main.o io.o fiber.o transport.o header.h
	nvcc main.o io.o fiber.o transport.o -o MCML_GPU
%.o: %.cu header.h
	nvcc -c -std=c++11  $<
clean:
	rm -rf main.o io.o fiber.o transport.o MCML_GPU

