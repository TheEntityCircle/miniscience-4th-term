GMSH_INCLUDE_PATH=/home/anganar/gmsh-4.9.4-Linux64-sdk/include/
GMSH_LIB_PATH=/home/anganar/gmsh-4.9.4-Linux64-sdk/lib/

all:
	for number in 1 2 3 4; do \
		g++ -I$(GMSH_INCLUDE_PATH) -L$(GMSH_LIB_PATH) sample0$$number.cpp -lgmsh -o sample0$$number; \
	done

clean:
	for number in 1 2 3 4; do \
        rm -f sample0$$number; \
    done
	rm -f *msh
